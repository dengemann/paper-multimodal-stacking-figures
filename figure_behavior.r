#'---
#'title: "3.1 Brain-age and behavior"
#'author: "Denis A. Engemann"
#'date: "9/1/2019"
#'output:
#'    html_document:
#'        code_folding:
#'            hide
#'    md_document:
#'        variant:
#'            markdown_github
#'---

#'In this series of visualizations we will explore how MRI, fMRI and #'MEG 
#'relate to behavioral data
#'
#'Questions we want to address here are
#'
#' 1. Are brain-age predictions associated with neuropsychological data?
#
#' 2. Do these associations go beyond the correlation of the neuropsychological
#'    data with age?
#'
#' 3. Do MEG, fMRI and MRI show similar relationships with behavior and
#'    does the strength of association change when data is combined.

#+ config
library(ggplot2)
library(ggbeeswarm)
library(ggrepel)
library(boot)
library(memoise)

# imports `color_cats` and `get_label_from_marker`
source('./utils.r')
source('./config.r')

PREDICTIONS <- './input_data/age_stacked_predictions_megglobal.csv'
BEHAVIOR <- './behavioral_data/data_summary.csv'

#' We will first read in the behavioral data.

#+ get_behavioral_data
data_3000 <- read.csv('behavioral_data/homeint_3000.tsv', sep = '\t')
data_add <- read.csv('behavioral_data/additional_3000.tsv', sep = '\t')
data_self <- read.csv(
    'behavioral_data/self_completion_questionnaire_3000.tsv', sep = '\t')
data_npsych <- read.csv("input_data/neuropsych_scores.csv")

names(data_npsych)[1] <- "Observations"

data_behavior <- Reduce(
    function(x, y) merge(x, y, by = "Observations", all.x = T),
    list(data_3000, data_add, data_self, data_npsych)
)

#' now the predictions

stacked_selection <- c(
  "ALL",
  "ALL_MRI",
  "ALL_no_fMRI",
  "MRI",
  "fMRI",
  "MEG_all"
)

data_pred_wide <- read.csv(PREDICTIONS)
data_pred_stacked_sel <- preprocess_prediction_data(
  df_wide = data_pred_wide, stack_sel = stacked_selection, drop_na = T)
names(data_pred_stacked_sel) <- sub(
    "X", "Observations", names(data_pred_stacked_sel))
names(data_pred_wide) <- sub(
    "X", "Observations", names(data_pred_wide))

data_behavior <- subset(data_behavior,
                        Observations %in% data_pred_wide$Observations)

age_data <- aggregate(age ~ Observations, data_pred_wide, unique)
data_behavior <- merge(data_behavior, age_data, by = "Observations")
stopifnot(sum(!data_behavior$age == age_data$age) == 0)


data_behavior <- within(data_behavior,
  {
    hours_in_bed[hours_in_bed > 20] <- NA
    hours_slept[hours_slept > 20] <- NA
    psqi[psqi > 20] <- NA
    HADS_anxiety[HADS_anxiety > 20] <- NA
    HADS_depression[HADS_depression > 20] <- NA 
    acer[acer > 500] <- NA
  }
)

# # this is what we may be interested in.
behavior_vars <- c(
 'psqi',  # not really continous (NRC), score
 'hours_slept', # count data
 'HADS_depression', # NRC, score
 'HADS_anxiety',  # NRC, score
 'acer', # NRC, score
 'mmse_i' # NRC, score
)

#' First thing to explore is a facet plot where each cell is one variable.
#' and scores are pitted against age predictions, grouped and coloured by
#' by brain age modality. The same thing is to be repeated after regressing out
#' age from any of these variables.
#'
#' In terms of preprocessing, we must compute the barin-age $\delta$
#' and also consider de-confounding the scores.
#'

data_pred_stacked_sel$delta <- with(data_pred_stacked_sel, pred - age)
data_pred_agg <- aggregate(delta ~ marker * Observations, data_pred_stacked_sel,
                           FUN = mean)

neuropsych_vars <- c(
  'BentonFaces',
  'CardioMeasures',
  'Cattell',
  'EkmanEmHex',
  'EmotionRegulation',
  'EmotionalMemory',
  'FamousFaces',
  'ForceMatching',
  'Hotel',
  'MotorLearning',
  'PicturePriming',
  'Proverbs',
  'RTchoice',
  'RTsimple',
  'Synsem',
  'TOT',
  'VSTMcolour')

selection <- Reduce(
  c,
  lapply(
    neuropsych_vars,
    function(x) {
      mask <- grepl(x, names(data_behavior))
      return(names(data_behavior)[mask])
}))

selection <- c("Observations", "age", selection, behavior_vars)

data_behavior <- data_behavior[,selection]
data_behavior_long <- reshape(
  data_behavior,
  varying = selection[-c(1, 2)],
  times = selection[-c(1, 2)],
  direction = "long",
  timevar = "psych_marker",
  v.names = "value"
)
data_behavior_long$psych_marker <- factor(data_behavior_long$psych_marker)
data_behavior_long$psych_family <- factor(
  sapply(strsplit(as.character(data_behavior_long$psych_marker), split = "_"),
         "[[", 1))

rename_marker <- function(x){
  x <- sapply(strsplit(x, split = "_"),
            abbreviate)
  x <- sapply(x, function(x) paste(x, collapse = ' '))
  return(x)
}

new_names <- setNames(
    rename_marker(levels(data_behavior_long$psych_marker)),
    levels(data_behavior_long$psych_marker))

my_labeller <- labeller(psych_marker = new_names)

fig3a_supp1 <- ggplot(
  data = data_behavior_long,
  mapping = aes(x = age, y = value, color = psych_marker)) +
  geom_point(alpha = 0.1, size = 0.5, show.legend = F) +
  geom_smooth(show.legend = F) +
  facet_wrap(~psych_marker, scales = "free",
             labeller = my_labeller) +
  mini_theme()
print(fig3a_supp1)

fname <- "./figures/elements_fig3_supplement1_age_functions."
ggsave(paste0(fname, "pdf"), plot = fig3a_supp1,
       width = save_width * 1.3,
       height = save_height * 1.3, useDingbats = F)
embedFonts(file = paste0(fname, "pdf"), outfile = paste0(fname, "pdf"))
ggsave(paste0(fname, "png"),
       width = save_width * 1.3,
       height = save_height * 1.3,
       plot = fig3a_supp1, dpi = 300)
knitr::include_graphics(paste0(fname, "png"), dpi = 200)

data_behavior_long_dec <- do.call(rbind, by(
  data_behavior_long,
  data_behavior_long$psych_marker,
  function(df){
    out <- data.frame(
      value = rep(NA, nrow(df)),
      Observations = df$Observations,
      psych_family = df$psych_family,
      psych_marker = df$psych_marker,
      age = df$age)
    good_mask <- !is.na(df$value)
    out$value[good_mask] <- resid(lm(value ~ poly(age, degree = 3),
                                  data = df[good_mask,]))
    return(out)
  }
))

fig3a_supp2 <- ggplot(
  data = data_behavior_long_dec,
  mapping = aes(x = age, y = value, color = psych_marker)) +
  geom_point(alpha = 0.1, size = 0.5, show.legend = F) +
  geom_smooth(show.legend = F) +
  facet_wrap(~psych_marker, scales = "free",
             labeller = my_labeller) +
  mini_theme()
print(fig3a_supp2)

fname <- "./figures/elements_fig3_supplement2_age_functions_deconfounded."
ggsave(paste0(fname, "pdf"), plot = fig3a_supp2,
       width = save_width * 1.3,
       height = save_height * 1.3, useDingbats = F)
embedFonts(file = paste0(fname, "pdf"), outfile = paste0(fname, "pdf"))
ggsave(paste0(fname, "png"),
       width = save_width * 1.3,
       height = save_height * 1.3,
       plot = fig3a_supp2, dpi = 300)
knitr::include_graphics(paste0(fname, "png"), dpi = 200)


data_full <- merge(data_pred_agg, data_behavior_long_dec, by = "Observations")

bootstrap_brain_age <- function(df) {
  set.seed(42)

  boot_result <- boot(
      within(df, {
    value <- scale(value)
    delta <- scale(delta)
  }),
      function(data, idx) {
        coef(lm(value ~ delta, data[idx,]))[['delta']]
      },
      R = 2000,
      parallel = "multicore",
      ncpus = 4)

  marker <- as.character(unique(df$marker)[1])
  psych_marker <- as.character(unique(df$psych_marker)[1])
  cat(marker, psych_marker, "\n")
  out <- data.frame(
        marker = marker,
        psych_marker = psych_marker,
        theta_boot = boot_result$t)
  out$theta <- boot_result$t0
  return(out)
}

memo_bootstrap_brain_age <- memoise(bootstrap_brain_age)

brain_age_boot <- do.call(rbind, by(
  data_full,
  list(data_full$psych_marker,
       data_full$marker),
  memo_bootstrap_brain_age
))

old_levels_ <- levels(brain_age_boot$marker)
levels(brain_age_boot$marker) <- c(
  "MRI, fMRI, MEG", "MRI, fMRI", "MRI, MEG", "fMRI", "MEG", "MRI")

colors_fig3a <- setNames(
  with(color_cats, c(black, orange, `blueish green`, blue, violet, vermillon)),
  c('MRI, fMRI, MEG', 'MRI, fMRI', 'MRI, MEG', 'MRI', 'fMRI', 'MEG'))


levels(brain_age_boot$psych_marker) <- new_names
fig3a <- ggplot(
  # data = subset(brain_age_boot, marker == "ALL"),
  data = brain_age_boot,
  mapping = aes(x = reorder(psych_marker, theta_boot, FUN = mean),
                y = theta_boot,
                fill = marker,
                color = marker)) +
  coord_flip() +
  stat_summary(geom = "boxplot", fun.data = my_quantiles, show.legend = F,
               alpha = 0.5, width = 0.8) +
  stat_summary(geom = "errorbar", fun.data = my_quantiles, show.legend = F) +
  geom_hline(yintercept = 0, color = 'black', linetype = 'dashed') +
  ylab(expression(beta[boot])) +
  xlab("Neuropsychological assessment") +
  theme(axis.text = element_text(size = 10)) +
  facet_wrap(~marker, nrow = 1) +
  scale_color_manual(values = colors_fig3a) +
  scale_fill_manual(values = colors_fig3a)
print(fig3a)

fname <- "./figures/elements_fig3a."
ggsave(paste0(fname, "pdf"), plot = fig3a,
       width = save_width * 1.5,
       height = save_height * 1.5, useDingbats = F)
embedFonts(file = paste0(fname, "pdf"), outfile = paste0(fname, "pdf"))
ggsave(paste0(fname, "png"),
       width = save_width * 1.5,
       height = save_height * 1.5,
       plot = fig3a, dpi = 300)
knitr::include_graphics(paste0(fname, "png"), dpi = 200)


data_full$marker <- factor(data_full$marker)
brain_age_pval <- do.call(rbind, by(
  data_full,
  list(data_full$psych_marker,
       data_full$marker),
  function(df) {
    df <- within(df, {
      value <- scale(value)
      delta <- scale(delta)
    })
    fit <- lm(value ~ delta, df)
    out <- data.frame(
      marker = unique(df$marker)[1],
      psych_marker = unique(df$psych_marker)[1],
      pval = summary(fit)$coefficients[2, 4],
      beta = coef(fit)[['delta']])
    return(out)
  }
))

levels(brain_age_pval$marker) <- levels(brain_age_boot$marker)
new_names2 <- setNames(
    LABELS_NPSYCH,
    levels(data_behavior_long$psych_marker))
levels(brain_age_pval$psych_marker) <- new_names2

percent_rank <- function(x) trunc(rank(x)) / length(x)
data_top_pval <- subset(
  brain_age_pval,
  pval <= 0.05
# -log10(brain_age_pval$pval) > -log10(0.05)
)

brain_age_pval <- within(brain_age_pval, {
  p_level <- rep('1', length(pval))
  p_level[pval <= 0.05] <- '2'
  # p_level[pval <= 0.01] <- 3
  # p_level[pval <= 0.001] <- 4
  p_level <- as.factor(p_level)
})

data_top_pval2 <- subset(
  brain_age_pval,
# percent_rank(-log10(brain_age_pval$pval)) >= .9
  p_level == 2
# -log10(brain_age_pval$pval) > -log10(0.05)
)

pos <- position_jitterdodge(jitter.width = 1.7, seed = 42,
                            dodge.width = 0.2)

fig3b <- ggplot(
  data = brain_age_pval,
  mapping = aes(y = -log10(pval),
                x = reorder(marker, -log10(pval), FUN = max),
                color = marker,
                shape = p_level,
                alpha = -log(pval),
                size = -log10(pval))) +
  geom_point(position = pos, show.legend = F) +
  # geom_label_repel(
  geom_text_repel(
    data = data_top_pval2,
    aes(label = psych_marker, size = -log10(pval)),
    parse = T,
    force = 2,
    position = pos,
    show.legend = F,
    segment.alpha = 0.4,
    size = 4
    ) +
  scale_alpha_continuous(range = c(0.05, 0.9), trans = 'sqrt', guide = F) +
  scale_shape_discrete(guide = F) +
  scale_color_manual(
    label = names(colors_fig3a),
    breaks = names(colors_fig3a),
    values = colors_fig3a) +
  scale_size_continuous(guide = FALSE) +
  guides(alpha = element_blank(), size = element_blank()) +
  geom_hline(yintercept = -log10(0.05)) +
  geom_hline(yintercept = -log10(0.01), linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.001), linetype = 'dotted') +
  # scale_x_discrete(breaks = NULL)+
  coord_cartesian(clip = 'off', ylim = c(-log10(1), 5)) +
  xlab(element_blank()) +
  ylab(expression(-log[10](p))) +
  guides(color = guide_legend(
          title = element_blank(),
          position = "top", nrow = 2,
           title.position = "left")) +
  theme(
      # legend.position = c(0.01, 0.9),
      axis.text = element_text(size = 14),
      legend.justification = 'left',
      legend.box = "horizontal")
print(fig3b)

fname <- "./figures/elements_fig3b."
ggsave(paste0(fname, "pdf"), plot = fig3b,
       width = save_width * 1.2,
       height = save_height, useDingbats = F)
embedFonts(file = paste0(fname, "pdf"), outfile = paste0(fname, "pdf"))
ggsave(paste0(fname, "png"),
       width = save_width * 1.2,
       height = save_height,
       plot = fig3b, dpi = 300)
knitr::include_graphics(paste0(fname, "png"), dpi = 200) 

#' Idea: do horizontal coef plot this way. Use p-values for ggploting to
#' indicate direction of effect. Change marker size /shape by sig level.
pos <- position_jitterdodge(jitter.width = 1.7, 
                            dodge.width = 0.2, 
                            seed = 42)
fig3c <- ggplot(
  data = brain_age_pval,
  mapping = aes(y = beta,
                x = reorder(marker, - log10(pval), FUN = max),
                color = marker,
                group = marker,
                alpha = -log(pval),
                size = -log(pval),
                shape = p_level)) +
  ylim(-.2, .2) +
  geom_point(position = pos, show.legend = F) +
  scale_color_manual(
    label = names(colors_fig3a),
    breaks = names(colors_fig3a),
    values = colors_fig3a) +
  ylab(expression(beta)) +
  xlab(element_blank()) +
  scale_size_continuous(guide = FALSE) +
  scale_alpha_continuous(range = c(0.05, 0.9), trans = 'sqrt', guide = F) +
  geom_text_repel(
    data = data_top_pval2,
    aes(label = psych_marker, size = -log10(pval)),
    parse = T,
    force = 2,
    position = pos,
    show.legend = F,
    segment.alpha = 0.4,
    size = 4
    ) +
  theme(
      axis.text = element_text(size = 14))

print(fig3c)

fname <- "./figures/elements_fig3c."
ggsave(paste0(fname, "pdf"),
       plot = fig3c,
       width = save_width * 1.2,
       height = save_height, useDingbats = F)
embedFonts(file = paste0(fname, "pdf"), outfile = paste0(fname, "pdf"))
ggsave(paste0(fname, "png"),
       width = save_width * 1.2,
       height = save_height,
       plot = fig3c, dpi = 300)

knitr::include_graphics(paste0(fname, "png"), dpi = 200) 

#+ session_info
print(sessionInfo())