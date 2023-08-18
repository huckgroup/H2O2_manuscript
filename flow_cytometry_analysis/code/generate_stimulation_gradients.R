###############################################
# profile object contains
# N - timepoints, N-1 - transitions
# timepoint: time, concentration, reactor_volume
# transition: added stimulus_volume, sample_volume, concentration_change
# attributes: ...

# procedure:
#   1) calculate exact profile
#   2) approximate profile with discrete flowrate values
#   3) calculate differences


#library(tidyverse)
library(purrr)
library(readr)
library(stringr)
library(ggplot2)

OUTPUT_ENCODING <- "latin1" # "ISO-8859-1" # "" #

#' f-n calculates flowrates for addition of stimulant and sampling to create gradient specified by `concentration_fn`
#' @param concentration_fn - a single parameter function
#' @param t_max - an end time of gradient; start time is always 0
#' @param t_unit - time units
#' @param V_start - reactors volume at the beginning (V(t = 0))
#' @param V_unit - volume units
#' @param df_sampling - data frame describing sampling; it consists of 3 vectors - time, volume, rate; @see <generation functions> for generation of sampling data frame
#' @param c_stimulus - stimulus concentration in a solution being added to a reactor
#' @returns gradient - data frame containing ....

gradient <- function(concentration_fn,
                     t_max,
                     t_unit = "s", #"s, ms, min"
                     V_start,
                     V_unit = "L",
                     df_sampling,
                     c_stimulus){
    if (t_unit != "s") {
        conv <- get_time_scaling_coef(t_unit)
        gradient(function(t) { concentration_fn(t / conv) },
                 t_max = t_max * conv,
                 t_unit = "s",
                 V_start,
                 V_unit,
                 df_sampling,
                 c_stimulus)
    } else if (V_unit != "L") {
        conv <- get_volume_scaling_coef(V_unit)
        gradient(concentration_fn,
                 t_max = t_max,
                 t_unit = t_unit,
                 V_start = V_start * conv,
                 V_unit = "L",
                 df_sampling,
                 c_stimulus = c_stimulus / conv)
    } else {
        gradient4(concentration_fn,
                  t_max,
                  dt = 0.1,
                  V_start,
                  df_sampling,
                  c_stimulus)
    }
}

gradient4 <- function(concentration,
                      tmax,
                      dt = 0.1,
                      Vstart,
                      df_sampling,
                      c_stimulus) {

    time <- seq(0,tmax, by = dt)
    c <- concentration(time)
    N <- length(c)
    dS <- vector(mode = "numeric", length = N - 1)
    dV <- vector(mode = "numeric", length = N - 1)
    V  <- vector(mode = "numeric", length = N)
    step_index <-
             vector(mode = "numeric", length = N - 1)
    V[1] <- Vstart
    index <- 1

    # the loop below iterates to solve for dV based on eq 9
    while(index < N) {
        c_now <- c[[index]]
        t_now <- time[[index]]
        v_now <- V[[index]]
        dc <- c[[index + 1]] - c_now
        dS[[index]] <- sampling_rate(sum(dS), t_now, df_sampling) * dt
        step_index[[index]] <- sum(get_sampling_steps(t_now, df_sampling)) * (dS[[index]] > 0)
#        dV[[index]] <- (dc * (v_now - dS[[index]]) + c_now * dS[[index]]) /
#                                    (c_stimulus - dc)
        dV[[index]] <- (dc * (v_now - dS[[index]])) / (c_stimulus - c[[index + 1]]) #(c_stimulus - c_now - dc)
        V[[index + 1]] <- V[[index]] + dV[[index]] - dS[[index]]
        index <- index + 1
    }
    list(time          = time,
         concentration = c,
         volume        = V,
         rate          = dV / dt,
         sampling_rate = dS / dt,
         sampling_index = step_index,
         dV            = dV,
         dS            = dS,
         dt            = dt,
         start_volume  = Vstart,
         df_sampling   = df_sampling,
         c_stimulus    = c_stimulus,
         concentration_fn = concentration) %>%
    list(continuous = .)
}
get_time_scaling_coef <- function(unit) {
    if (unit == "s") {
        1
    } else if (unit == "min") {
        60
    } else if (unit == "h") {
        3600
    } else if (unit == "ms") {
        0.001
    } else {
        unit %>%
            str_c("unit `", ., "` not recognized, assuming values are in `s`", sep = "") %>%
            print
        1
    }
 }
get_rate_scaling_coef <- function(unit) {
    if (is.null(unit) || length(unit) == 0) {
        print("get_rate_scaling_coef(): NULL of zero length unit parameter, assuming L/s")
        1
    } else if (length(unit) > 1) {
        map_dbl(unit, get_rate_scaling_coef)
    } else {
        components <- str_split(unit, "/")[[1]]
        if (length(components) == 2) {
            get_volume_scaling_coef(components[[1]]) / get_time_scaling_coef(components[[2]])
        } else {
            print("get_rate_scaling_coef(): failed parsing rate unit, assuming L/s")
            1
        }
    }
}
get_volume_scaling_coef <- function(unit) {
    unit_ <- gsub("l", "L", unit)
    if (unit_ == "L") {
        1
    } else if (unit_ == "mL") {
        0.001
    } else if (unit_ == "mkL" || unit_ == "μL") {
        0.000001
    } else {
        unit %>%
            str_c("unit `", ., "` not recognized, assuming values are in `L`", sep = "") %>%
            print
        1
    }
}
get_sampling_steps <- function(t, df_sampling) {
    # f-n gets vector of sampling steps used so far (TRUE - used, FALSE - not used)
    # df_sampling - data frame with sampling infor
    #               [["time"]]   sampling time
    #               [["volume"]] sampling volume
    #               [["rate"]]   sampling rate
    # V            - total sampled volume up to time t
    # t            - current time
    if (nrow(df_sampling) < 1) {
        return(0) # no steps to sample: empty vector
    } else {
        df_sampling[["time"]] <= t
    }
}
sampling_rate <- function(V, t, df_sampling) {
    # f-n estimates if sampling has to be done (rate | 0)
    # df_sampling - data frame with sampling infor
    #               [["time"]]   sampling time
    #               [["volume"]] sampling volume
    #               [["rate"]]   sampling rate
    # V            - total sampled volume up to time t
    # t            - current time
    points_to_use <- get_sampling_steps(t,df_sampling)
    index <- sum(points_to_use)
    if (index < 1) {
        0
    } else if (V < sum(df_sampling[["volume"]][points_to_use])) {
#debug#        print("getting non-zerro rate")
        df_sampling[["rate"]][[index]]
    } else {
        0
    }
}
plot_gradient <- function(gradient, parameter = "rate", V_unit = "L", t_unit = "s") {
  value1 <- gradient[["continuous"]][[parameter]] * 1/get_volume_scaling_coef(V_unit)
  time1 <- gradient[["continuous"]][["time"]][seq(1, length(value1))] / get_time_scaling_coef(t_unit)
  label1 <- rep(parameter, length(value1))

  value2 <- gradient[["approximation"]][[parameter]]
  time2 <- gradient[["approximation"]][["time"]][seq(1, length(value2))]
  label2 <- rep(paste0(parameter, "_approx"), length(value2))

  data.frame(x = c(time1, time2),
             y = c(value1, value2),
             lbl = c(label1, label2)) %>%
  ggplot(.)+ #replace 1 to gradients comp name!!!
    geom_line(mapping = aes(x = x,
                            y = y,
                            color = lbl)) +
    ggtitle(paste0("custom gradient: ", parameter)) +
    xlab("time") +
    ylab(parameter)
}
empty_df_sampling <- function() {
    tibble(time = list(), volume = list(), rate = list())
}
# --- create df_sampling from input string ---
#       for a single sampling point end the string with `row_sep`
#       (e.g. use "100,0.1,10000\n" instead of "100,0.1,10000",
#       "\n" is row_sep in this example)
create_df_sampling <- function(input, row_sep = "\n", col_sep = ",", units = c("s", "L", "uL/h")) {
    if (is.null(input) ||
        length(input) == 0 ||
        (length(input) == 1 && input == "")) {
        empty_df_sampling()
    } else if (row_sep != "\n") {
        create_df_sampling(gsub(row_sep, "\n", input), col_sep = col_sep)
    } else if (col_sep != ",") {
        create_df_sampling(gsub(col_sep, ",", input), row_sep = row_sep)
    } else if (length(input) > 1) {
        print("multiple inputs")
        print(input)
        map(input, ~create_df_sampling(.x, row_sep = row_sep, col_sep = col_sep))
    } else {
        df <- read_csv(input, col_names = c("time", "volume", "rate"))
        convert_df_scaling_units(df, units)
    }
}
create_df_sampling_low <- function(time, volume, rate, units) {
    if (any(map_lgl(list(time, volume, rate, units), is.null)) ||
        !(length(time) == length(volume) &&
          (length(time) == length(rate) || length(rate) == 1) &&
          length(units) == 3)) {
        print("create_df_sampling_low(): incorrect input format; empty df_sampling created")
        empty_df_sampling()
    } else {
        rate_ <- if (length(rate) == 1) rep(rate, length(time)) else rate
        data.frame(time = time, volume = volume, rate = rate_) %>%
            convert_df_scaling_units(., units)
    }
}
create_df_sampling_uniform <- function(N, d_sample, dt, rate, units) {
    if (N < 1) {
        empty_df_sampling()
    } else {
        create_df_sampling_low(time   = seq(1, N) * dt,
                               volume = rep(d_sample, N),
                               rate   = rep(rate, N),
                               units  = units)
    }
}
convert_df_scaling_units <- function(df_sampling_raw, units) {
    conv_t <- get_time_scaling_coef(units[[1]])
    conv_V <- get_volume_scaling_coef(units[[2]])
    conv_rate <- get_rate_scaling_coef(units[[3]])

    print(df_sampling_raw)
    print(conv_t)
    print(conv_V)
    print(conv_rate)

    data.frame(time   = df_sampling_raw[["time"]] * conv_t,
               volume = df_sampling_raw[["volume"]] * conv_V,
               rate   = df_sampling_raw[["rate"]] * conv_rate)
}
assign_sampling_index <- function(x, N) {
    map_dbl(x,
            function(x_val) {
                if (x_val > 0) {
                    (x_val - 1) %% N + 1
                } else {
                    0
                }
            })
}
create_sampling_rates <- function(gradient, n_sampling_pumps = 2) {
    g <- gradient[["continuous"]]
    sampling_rate_general <- g[["sampling_rate"]]
    sampling_index <- assign_sampling_index(g[["sampling_index"]], n_sampling_pumps)

    # print(sampling_index)
    # plot(sampling_index)
    # print(table(sampling_index))
    max_val <- max(sampling_index)
    print(str_c("max_val: ", max_val, sep = ""))
    if ( max_val > 0) {
        # constant `-1` ensures negative flow rates (loading into syringe)
        map(seq(1,max_val), ~ -1 * sampling_rate_general * (sampling_index == .x))
    } else {
        list()
    }
}
export_nemesys_profile_step_v4 <- function(file, rate_vector, dt, coef_time = 1000, coef_rate = 3.6e9, units = "µl/h", index = 0, digits = 2, decimal.mark = ",") {
# f-n exports calculated gradient
#       assumes that calculation are done in liters (L)  and seconds (s),
#       default coeficients convert rates to uL/h and intervals to ms
# µl/h
  f <- file(file,
            open = "wb", #if (append) "a" else "w",
            encoding = OUTPUT_ENCODING)
#  open(f)

  my_write <- function (x, append = TRUE) {
    cat(x,
        file = f,
        #append = append,
        sep = "\n")
  }
#debug#  print(units)
  my_write(units, FALSE)
#debug#  print(index)
  my_write(index)

  rate_duration <-
      rate_vector %>%
      `*`(coef_rate) %>%
      round(digits = digits) %>%
      rle

    print(paste(c("Amount (\u03BCL) of stimulus needed: ", round(sum(rate_duration[[1]]/36000*rate_duration[[2]]), digits = 2)), collapse = ""))


  extra_step <- data.frame(time = 100 * coef_time,
                           rate = 0,
                           valve = 255)

  data.frame(time = (rate_duration[["lengths"]] * dt * coef_time),
             rate = rate_duration[["values"]],
             valve = rep(255, length(rate_duration[["lengths"]]))) %>%
      rbind(extra_step) %>%  # add extra step for pause at the end
                             # to prevent immediate recycling of a
                             # profile step
      format(scientific = F, decimal.mark = decimal.mark) %>%
      #gsub(" ", "", .)       %>%
      write.table(file = f,
                  sep = "\t",
                  #append = TRUE,
                  col.names = FALSE,
                  row.names = FALSE,
                  #dec = ",",
                  eol = "\n",
                  quote = FALSE,
                  fileEncoding = OUTPUT_ENCODING)
  close(f)
}

export_nemesys_profile_step_v4_experimental <- function(file, rate_vector, dt, coef_time = 1000, coef_rate = 3.6e9, units = "µl/h", index = 0, digits = 2) {
# f-n exports calculated gradient
#       assumes that calculation are done in liters (L)  and seconds (s),
#       default coeficients convert rates to uL/h and intervals to ms
# µl/h
  f <- file(file,
            open = "wb", #if (append) "a" else "w",
            encoding = OUTPUT_ENCODING)
#  open(f)

  my_write <- function (x, append = TRUE) {
    cat(x,
        file = f,
        #append = append,
        sep = "\n")
  }
#debug#  print(units)
  my_write(units, FALSE)
#debug#  print(index)
  my_write(index)

  rate_duration <- rle(get_approximated_rate(rate_vector, coef_rate, digits))
#      rate_vector %>%
#      `*`(coef_rate) %>%
#      round(digits = digits) %>%
#      rle

  data.frame(time = (rate_duration[["lengths"]] * dt * coef_time),
             rate = rate_duration[["values"]],
             valve = rep(255, length(rate_duration[["lengths"]]))) %>%
      format(scientific = F, decimal.mark = ",") %>%
      #gsub(" ", "", .)       %>%
      my_write(.)

  close(f)
}
export_nemesys_profile_gradient <- function(file,
                                            gradient,
                                            coef_time = 1000,
                                            coef_rate = 3.6e9,
                                            units = "μl/h",
                                            index = 0,
                                            digits = 2,
                                            n_sampling_pumps = 2,
                                            decimal.mark = ",") {
      make_export_path <- function(lbl) {
          folder <- dirname(file)
          b_name <- basename(file)
          new_name <- str_c(b_name, "_", lbl, ".nfp", sep = "")
          file.path(folder, new_name)
      }

      data_to_export <- reduce(list(create_sampling_rates(gradient, n_sampling_pumps)),
                               .init = list(gradient[["continuous"]][["rate"]]),
                               .f = append)
      # add extra step that does nothing for 100 seconds after profile is done
      # profiles are reused constantly and have to be stoped manually,
      # an extra step is to give time to user to stop experiment.

      data_to_export %>%

          str

      lbls <- "stimulus" %>%
                c(if (length(data_to_export) <= 1) {
                      NULL
                  } else {
                      seq(1, length(data_to_export) - 1)
                  })
      export_path <- map_chr(lbls, make_export_path)

#debug#      print("about to walk")
    walk2(export_path,
          data_to_export,
          ~export_nemesys_profile_step_v4(
               file = .x,
               rate_vector = .y,
               dt = gradient[["continuous"]][["dt"]],
               coef_time = coef_time,
               coef_rate = coef_rate,
               units = units,
               index = index,
               digits = digits,
               decimal.mark = decimal.mark))
}
calc_gradient <- function(v0, dt, rate_stimulant, rate_sampling, c_stimulant_0) {
    N <- length(rate_stimulant)
    conc <- rep(0, N + 1)
    volume <- rep(v0, N + 1)
    time <- seq(0, N) * dt
    for (index in seq(1,N)) {
        V_stimulant_added <- rate_stimulant[[index]] * dt
        V_sample  <- rate_sampling[[index]] * dt
        volume[[index + 1]] <- volume[[index]] + V_stimulant_added - V_sample
        n_stimulant <- volume[[index]] * conc[[index]] +
                       V_stimulant_added * c_stimulant_0 -
                       V_sample * conc[[index]]
        conc[[index + 1]] <- n_stimulant / volume[[index + 1]]
    }
    list(time = time, volume = volume, concentration = conc)
}
check_gradient <- function(gradient, coef_rate, digits) {
    g  <-  gradient[["continuous"]]
    df1 <- data.frame(time = g[["time"]],
                      concentration = g[["concentration"]],
                      volume = g[["volume"]],
                      type = rep("formula", length(g[["time"]])))


    V_init  <- g[["start_volume"]]
    dt  <- g[["dt"]]
    c_stim_0  <- g[["c_stimulus"]]
    check1 <- calc_gradient(V_init,
                            dt,
                            g[["rate"]],
                            g[["sampling_rate"]],
                            c_stim_0)
    check1[["type"]] <- rep("check_exact", length(check1[["time"]]))
    df2 <- as.data.frame(check1)

    check2 <- calc_gradient(V_init,
                            dt,
                            get_approximated_rate(g[["rate"]],
                                                  coef_rate,
                                                  digits) / coef_rate,
                            get_approximated_rate(g[["sampling_rate"]],
                                                  coef_rate,
                                                  digits) / coef_rate,
                            c_stim_0)
    check2[["type"]]  <- rep("approximation", length(check2[["time"]]))
    df3 <- as.data.frame(check2)

    rbind(df1, df2) %>%
        rbind(.,df3)
}

plot_gradient_check  <- function(gradient, coef_rate = 3.6e9, digits) {
    check_gradient(gradient, coef_rate, digits) %>%
        ggplot(aes(x = time, y = concentration, color = type)) +
        geom_line()
}

get_approximated_rate <- function(rate_vector, coef_rate, digits) {
      round(rate_vector * coef_rate, digits = digits)
}


