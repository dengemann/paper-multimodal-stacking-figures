#'---
#'title: "1.1 MEG performance"
#'author: "Denis A. Engemann"
#'date: "8/3/2019"
#'output:
#'    html_document:
#'        code_folding:
#'            hide
#'    md_document:
#'        variant:
#'            markdown_github
#'---

#'In this series of visualizations we will explore the linear AGE fits
#'and the stacked predictions for different MEG variables.
#'We shall answer the following questions: 
#'
#'1. What is the impact of stacking?
#'
#'2. Which type of MEG feature is most influential on the results?
#'
#'3. Is the feature-type or the frequency band more impotant?
#'
#'4. Does the stacker learn additive or interactive effects between features?

#+
library(ggplot2)
library(ggbeeswarm)
library(ggrepel)

# imports `color_cats` and `get_label_from_marker`
source('./utils.r')
source('./config.r')

#' We'll first set up the data and configure plotting.
#' Also, in this analysis we'll focus on the scores, not the predictions.

#+ setup
PREDICTIONS <- './input_data/age_stacked_predictions_megglobal.csv'
SCORES <- './input_data/age_stacked_scores_megglobal.csv'
IMPORTANCES <- './outputs/age_stacked_importance_8.csv'

data_pred_wide <- read.csv(PREDICTIONS)
data_scores_wide <- read.csv(SCORES)

names(data_scores_wide) <- fix_column_names(names(data_scores_wide))
names(data_pred_wide) <- fix_column_names(names(data_pred_wide))

# Move to long format.
reshape_names <- names(data_scores_wide)[!names(data_scores_wide) %in%
                                         c('repeat_', 'repeat_idx')]
data_scores <- reshape(data = data_scores_wide,
                       direction = 'long',
                       varying = reshape_names,  # values
                       v.names = 'MAE',
                       times = factor(reshape_names), # keys
                       timevar = 'marker',
                       new.row.names = NULL)

data_scores['is_meg'] <- grepl('MEG', data_scores$marker)
#'Remap the names to be more compact for plotting.

#+ stack_labels
stacked_keys <- c(
  "MEG_handcrafted",
  "MEG_powers",
  "MEG_powers_cross_powers",
  "MEG_powers_cross_powers_handrafted",
  "MEG_cat_powers_cross_powers_correlation",
  "MEG_cat_powers_cross_powers_correlation_handcrafted",
  "MEG_cross_powers_correlation",
  "MEG_powers_cross_powers_correlation",
  "MEG_all"
)

stacked_labels <- c(
    'O',
    'P',
    'P~XP[f]',
    'P~XP[f]~O',
    'P[cat]~XP[f]~C[f]',
    'P[cat]~XP[f]~C[f]~O',
    'XP[f]~C[f]',
    'P~XP[f]~C[f]',
    'P~XP[f]~C[f]~O'
)

#+ subset_meg
data_scores$prediction <- factor(ifelse(
  data_scores$marker %in% stacked_keys, 'stacked', 'linear'))

# For now ignore MRI
data_scores <- subset(data_scores, is_meg == T)
data_scores$marker <- factor(data_scores$marker)
data_scores$repeat_idx <- factor(data_scores$repeat_idx)

#'Now we compute best results, we need to aggregate by fold idx.

#+ agg_by_fold
data_scores_cv <- by(data_scores,
                     list(data_scores$repeat_idx,
                          data_scores$marker),
                     function(x){
                       data.frame(
                         marker = unique(x$marker),
                         repeat_idx = unique(x$repeat_idx),
                         cv_mean = mean(x$MAE),
                         cv_std = sd(x$MAE))
                     })
data_scores_cv <- do.call(rbind, data_scores_cv)

data_scores_cv$prediction <- factor(ifelse(
  data_scores_cv$marker %in% stacked_keys, 'stacked', 'linear'))

best_cv_mean_linear <- min(subset(data_scores_cv,
                                  prediction == 'linear')$cv_mean)
best_cv_mean_stacked <- min(subset(data_scores_cv,
                                   prediction == 'stacked')$cv_mean)

#'We are ready to plot the error distributions.

#+ fig1a
n_models <- by(data_scores,
               data_scores$prediction,
               function(x) data.frame(
                 marker = unique(x$prediction),
                 n_models =  length(unique(x$marker))))
n_models <- do.call(rbind, n_models)

labels_fig1 <- c(
  sprintf("%s (n=%d)", n_models[2,1], n_models[2,]$n_models),
  sprintf("%s (n=%d)", n_models[1, 1], n_models[1,]$n_models)
)

labels_fig1 <- setNames(labels_fig1, c('stacked', 'linear'))
colors_fig1 <- setNames(with(color_cats, c(black, `blueish green`)),
                        c('stacked', 'linear'))
data_scores$prediction <- relevel(data_scores$prediction, ref = "stacked")

fig1a <- ggplot(data = data_scores,
                mapping = aes(x = MAE,
                              color = prediction, fill = prediction)) +
  stat_density(trim = T, geom = 'area', size = 1, alpha = 0.3) +
  scale_x_log10(breaks = seq(0, 60, 5)) +
  scale_color_manual(
    values = colors_fig1,
    labels = labels_fig1) + 
  scale_fill_manual(
    values = colors_fig1,
    labels = labels_fig1) + 
  guides(
    color = guide_legend(nrow = 1, title.position = "left")) +
  theme(legend.position = c(0.4, 0.99)) +
  xlab("MAE (years)") +
  ylab("Density") +
  annotate(geom = "text", x = DUMMY_ERROR - 0.5, y = 5,
           label =  "predicting~bar(y)",
           parse = T,
           angle =90, size = annotate_text_size) +
  geom_vline(xintercept=DUMMY_ERROR, color = "black", linetype = 'dashed') +
  geom_rug(alpha = 0.5, mapping = aes(color = prediction, x = cv_mean),
           data = data_scores_cv)
print(fig1a)


fname <- "./figures/elements_fig1_error_distro."
ggsave(paste0(fname, "pdf"),
       width = save_width, height = save_height, useDingbats = F)
embedFonts(file = paste0(fname, "pdf"), outfile = paste0(fname, "pdf"))
ggsave(paste0(fname, "png"),
       width = save_width, height = save_height,
       dpi = 300)
knitr::include_graphics(paste0(fname, "png"), dpi = 200)

#'Remap the names to be more compact for plotting.

#'We can see that, across models, stacking stablizes and improves the results.
#'
#'Let's now look at the most interesting stacked models.
#'
#'However, we should first add some more description to the data.

#+ subset_meg_stack
data_meg_stacked <- subset(data_scores, prediction == 'stacked')

data_meg_stacked$family <- 'Combined'
data_meg_stacked$family[
  data_meg_stacked$marker == 'MEG_powers'] <- 'Power'
data_meg_stacked$family[
  data_meg_stacked$marker == 'MEG_cross_powers_correlation'] <- 'Connectivity'
data_meg_stacked$family[
  data_meg_stacked$marker == 'MEG_handcrafted'] <- 'Other'
data_meg_stacked$family[
  data_meg_stacked$marker == 'MEG_all'] <- 'Full'

data_meg_stacked$family <- factor(data_meg_stacked$family)

data_meg_stacked$marker <- factor(
  data_meg_stacked$marker,
  stacked_keys,
  stacked_labels)

#+ fig1b
color_breaks <- c('Full', 'Combined', 'Power', 'Connectivity', 'Other')
color_values <-setNames(
  with(color_cats, c(black, orange, vermillon, blue, `sky blue`)),
  color_breaks
)
color_labels <- tolower(c(color_breaks[1:4], "sensor"))
sel_idx <- data_meg_stacked$marker %in% c(
  "O",
  "P",
  "XP[f]~C[f]",
  "P[cat]~XP[f]~C[f]",
  "P~XP[f]~C[f]~O"
)

if(FALSE){  # toggle to select everything
  sel_idx <- rep(T, nrow(data_meg_stacked))
}

data_meg_stacked_sel <- data_meg_stacked[sel_idx,]
data_meg_stacked_sel$marker <- factor(data_meg_stacked_sel$marker)
sort_idx <- order(by(data_meg_stacked_sel,
                     data_meg_stacked_sel$marker,
                     function(x) mean(x$MAE)))

fig1b <- ggplot(
  data = data_meg_stacked_sel,
  mapping = aes(y = MAE, x = reorder(marker, MAE, mean),
                color = family, fill = family)) +
  coord_flip() +
  geom_hline(yintercept=DUMMY_ERROR, color = "black", linetype = 'dashed') +
  geom_boxplot(show.legend = T, outlier.shape = NA, alpha = 0.5, size = 0.7) +
  stat_summary(geom = 'text',
               mapping = aes(label  = sprintf("%1.1f", ..y..)),
               fun.y= mean, size = 3.2, show.legend = FALSE,
               position = position_nudge(x=-0.49)) +
  geom_beeswarm(alpha=0.3, show.legend = F, size = 2.4) +
  guides(
    color = guide_legend(nrow = 3, title.position = "top")) +
  ylab("MAE (years)") +
  xlab("MEG stacking models") +
  annotate(geom = "text", y = DUMMY_ERROR - 0.3, x = 2.5,
           label =  "predicting~bar(y)",
           parse = T,
           angle = 90, size = annotate_text_size) +
  # scale_x_discrete(
  #   labels = parse(text = levels(data_meg_stacked$marker)[sort_idx])) +
  scale_color_manual(values = color_values,
                     labels = color_labels,
                     breaks = color_breaks, name = "stacking") +
  scale_fill_manual(values = color_values,
                    labels = color_labels,
                    breaks = color_breaks, name = "stacking") +
  guides(color = guide_legend(nrow = 1, title.position = 'top'),
         fill = guide_legend(nrow = 1)) +
  scale_y_continuous(breaks = seq(3, 20, 1)) +
  theme(
        axis.text.y = element_blank(),
        legend.position = 'top',
        legend.justification = 'left',
        legend.text.align = 0)
print(fig1b)

fname <- "./figures/elements_fig1_stacking_models."
ggsave(paste0(fname, "pdf"),
       width = save_width, height = save_height, useDingbats = F)
embedFonts(file = paste0(fname, "pdf"), outfile = paste0(fname, "pdf"))
ggsave(paste0(fname, "png"),
       width = save_width, height = save_height,
       dpi = 300)
knitr::include_graphics(paste0(fname, "png"), dpi = 200)

#'Now we investigate variable importance. We need to do some data massaging
#'in order to have everyhing in a neat structure.

#+ process_importance
stack_models <- c(
    'MEG-all-no-diag',
    'MEG all'
)

importance_methods <- c(
  'permutation',
  'rf_msqrt',
  'rf_m1',
  'et_m1'
)

# we'll expand the grid of options
importance_params <- expand.grid(stack_models, importance_methods)

data_importance <- read.csv(IMPORTANCES)

importance_results <- vector("list", nrow(importance_params))
for (ii in seq_along(importance_results)){
  this_importance <- subset(
    data_importance,
    stack_model == as.character(importance_params[ii, 1]) &
    mod_type == as.character(importance_params[ii, 2]))
  importance_results[[ii]] <- get_importances(this_importance, data_scores)
}
importance_results <- do.call(rbind, importance_results)
#'Let's add some discriptors to enrich out GG capabilites.
#'We'll also make some nicer labels for plotting.

#+ enrich_importance_info
importance_results$family <- get_family_from_marker(importance_results$marker)
importance_results$variant <- get_variant_from_marker(importance_results$marker)

importance_results$marker_label <- get_label_from_marker(
  importance_results$marker)

importance_results$stack_model <- sub(
  "MEG all", "MEG_all", importance_results$stack_model)
#' We are ready to plot our importance results.

#+ fig1c
color_breaks <- c('power', 'cross', 'corr', 'other')
color_labels <- c('power', 'cov', 'corr', 'sensor')
color_values <-setNames(
  with(color_cats, c(vermillon, orange, blue, `sky blue`)),
  color_breaks
)


shape_breaks <- c('power', 'envelope', '1/f', 'peak')
shape_values <- setNames(
  c(15, 17, 18, 16),
  shape_breaks
)

method <- 'permutation'
stack_model <- 'MEG_all'
sub_data <-  importance_results[
  importance_results$method == method &
  importance_results$stack_model == stack_model,]

sub_data$variant <- factor(sub_data$variant,
  levels = c("power", "envelope", "1/f", "base"),
  labels=c("power", "envelope", "1/f", "peak"))

# lets dump the most important markers
percent_rank <- function(x) trunc(rank(x)) / length(x)
write.csv(
  sub_data[order(percent_rank(sub_data$importance), decreasing  = T),][1:10,],
  file = paste0(paste('./outputs/importances', method, stack_model, sep = '_'),
                '.csv'))

fig1c <- ggplot(
  data = sub_data,
  mapping = aes(x = importance,
                y = cv_score,
                color = family,
                shape = variant)) +
  geom_point(size = 3., alpha = 0.8) +
  scale_shape_manual(
    values = shape_values,
    breaks =  shape_breaks,
    labels = shape_breaks
    ) +
  scale_color_manual(
    values = color_values,
    labels = color_labels,
    breaks = color_breaks) +
  guides(
    color = guide_legend(
      nrow = 4, title.position = "top", order = 1),
    shape = guide_legend(
      nrow = 4, title.position = "top", order = 2)) +
  theme(legend.position = c(0.55, 0.85),
        legend.justification = 'left',
        legend.box = "horizontal",
        legend.direction = 'horizontal') +
  scale_y_log10(breaks = seq(0, 50, 5)) +
  scale_x_log10() +
  ylab("MAE (years)") +
  xlab("Variable importance") + 
  geom_label_repel(
    data = subset(sub_data, percent_rank(sub_data$importance) >= .9),
    aes(label = marker_label),
    parse = T,
    force = 2,
    show.legend = F,
    segment.alpha = 0.4,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  )
print(fig1c)

fname <- paste0("./figures/elements_fig1c_importance",
                "_", stack_model, "_", method, ".")
ggsave(paste0(fname, "pdf"),
       width = save_width, height = save_height, useDingbats = F)
embedFonts(file = paste0(fname, "pdf"), outfile = paste0(fname, "pdf"))
ggsave(paste0(fname, "png"), width = save_width, height = save_height,
       dpi = 300)
knitr::include_graphics(paste0(fname, "png"), dpi = 200)

#'Time to take a look at the partial dependencies.
#'We'll also enrich the data with descriptors for plotting, as above.

#+ read_dependence
OUT_DEPENDENCE_1D <- './input_data/age_stacked_dependence_model-full-1d.csv'
data_dependence_1d <- read.csv(OUT_DEPENDENCE_1D)
data_dependence_1d$marker <- factor(data_dependence_1d$marker)
data_dependence_1d$model_type <- factor(data_dependence_1d$model_type)
data_dependence_1d$model <- factor(data_dependence_1d$model)

#+ fig1d
data_dependence_1d$marker <- fix_marker_names(data_dependence_1d$marker)
data_dependence_1d$family <- get_family_from_marker(data_dependence_1d$marker)
data_dependence_1d$variant <- get_variant_from_marker(data_dependence_1d$marker)
data_dependence_1d$marker_label <- get_label_from_marker(
  data_dependence_1d$marker)

color_breaks_1d <- c('power', 'cross', 'corr')
color_labels_1d <- c('power', 'cov', 'corr')
color_values_1d <-setNames(
  with(color_cats, c(vermillon, orange, blue)),
  color_breaks_1d
)

data_sel <- subset(data_dependence_1d,
  model_type == 'rf_msqrt' & model == 'MEG all')

data_sel_max <- aggregate(. ~ marker * family * variant * marker_label,
                          data = data_sel, FUN = max)

fig1d <- ggplot(data = data_sel,
       mapping = aes(x = value,  y = pred,
                     group = marker, color = family,
                     linetype = variant)) +
  geom_line(size = 1) + 
  xlab('Input age (years)') +
  ylab('Age prediciton (years)') +
  scale_color_manual(
    values = color_values_1d,
    labels = color_labels_1d,
    breaks = color_breaks_1d) +
  guides(
    color = guide_legend(
      nrow = 4, title.position = "top", order = 1),
    shape = guide_legend(
      nrow = 4, title.position = "top", order = 2))+
    theme(
      legend.position = c(0.01, 0.8),
        legend.justification = 'left',
        legend.box = "horizontal",
      ) +
  geom_label_repel(
    data = data_sel_max,
    aes(label = marker_label),
    parse = T,
    force = 2,
    show.legend = F,
    segment.alpha = 0.4,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  )
print(fig1d)  

fname <- "./figures/elements_fig1d_dependence."
ggsave(paste0(fname, "pdf"), width = save_width,
       height = save_height, useDingbats = F)
embedFonts(file = paste0(fname, "pdf"), outfile = paste0(fname, "pdf"))
ggsave(paste0(fname, "png"), width = save_width,
       height = save_height, dpi = 300)
knitr::include_graphics(paste0(fname, "png"), dpi = 200)


#'Let us also analyse the 2D partial dependence between the best 4 markers.

#+ pdp_2d
# is already in long format.
OUT_DEPENDENCE_2D <- './input_data/age_stacked_dependence_model-full-2d.csv'
data_dependence2d <- read.csv(OUT_DEPENDENCE_2D)
data_dependence2d$marker <- fix_marker_names(data_dependence2d$marker)
data_dependence2d$var_x <- fix_marker_names(data_dependence2d$var_x)
data_dependence2d$var_y <- fix_marker_names(data_dependence2d$var_y)
data_dependence2d$model <- sub(" ", "_", data_dependence2d$model)

# get label
data_dependence2d$x_label <- get_label_from_marker(
  data_dependence2d$var_x)
data_dependence2d$y_label <- get_label_from_marker(
  data_dependence2d$var_y)
# ... family
data_dependence2d$x_family <- get_family_from_marker(
  data_dependence2d$var_x)
data_dependence2d$y_family <- get_family_from_marker(
  data_dependence2d$var_y)
# ... variant
data_dependence2d$x_variant <- get_variant_from_marker(
  data_dependence2d$var_x)
data_dependence2d$y_variant <- get_variant_from_marker(
  data_dependence2d$var_y)

pdp2d_cases <- c(
  "envelope_diag--mne_envelope_cross_alpha",
  "envelope_diag--mne_envelope_cross_beta_low",
  "envelope_diag--power_diag",
  "mne_envelope_cross_alpha--mne_envelope_cross_beta_low",
  "power_diag--mne_envelope_cross_alpha",
  "power_diag--mne_envelope_cross_beta_low"
)

pdp2d_labels <- c(
  "E[cat]-alpha[env~x]",
  "E[cat]-beta[low~env~x]",
  "E[cat]-P[cat]",
  "alpha[env~x]-beta[low~env~x]",
  "P[cat]-alpha[env~x]",
  "P[cat]-beta[low~env~x]"
)

#  make nicer marker labels for titles.
data_dependence2d$marker_label <- data_dependence2d$marker
for (i_case in seq_along(pdp2d_cases)) {
  data_dependence2d$marker_label <- sub(
    pdp2d_cases[[i_case]], pdp2d_labels[[i_case]],
    data_dependence2d$marker_label)
}

data_dependence2d_sel <- subset(
  data_dependence2d,
  model == "MEG_all" & model_type == 'rf_msqrt'
)

this_data <- data_dependence2d[
  data_dependence2d$model == "MEG_all" &
  data_dependence2d$model_type == 'rf_msqrt',]

this_breaks <- seq(49, 59, 1)

fig1e <- ggplot(data = this_data,
                mapping = aes(x = x, y = y, z = pred)) +
    geom_raster(aes(fill = pred), show.legend = T) +
    stat_contour(breaks = this_breaks,
                 color = "white", bins = 7, show.legend = F) +
    scale_fill_viridis_c(breaks = this_breaks,
                         name = expression(hat(y)),
                         guide = guide_colorbar(barheight = 10)) +
    theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
    facet_wrap(~marker_label, scales = "free", labeller = label_parsed) +
    xlab("Input age 1") +
    ylab("Input age 2")

fname <- paste0("./figures/elements_fig1e_meg_pdp_2d.")
ggsave(paste0(fname, "pdf"), plot = fig1e,
      width = save_width, height = save_height, useDingbats = F)
embedFonts(file = paste0(fname, "pdf"), outfile = paste0(fname, "pdf"))
ggsave(paste0(fname, "png"), plot = fig1e,
        width = save_width, height = save_height,
      dpi = 300)
knitr::include_graphics(paste0(fname, "png"), dpi = 200)

#'We can see at least two things. 1) The random forest combines the age input
#'rather additively, as can be seen by the diagonalish pattern when moving
#'through age input ranges. 2)  Howver, the patterns is  not always perfectly
#'diagonal, which implies that some of  the inputs are kicking later or earlier.
#'There is also a somewhat step-like transition, which makes somewhat
#'rectangular surfaces. In total, this makes additive but non-linear
#'integration.

fname_fig1 <- './figures/fig_1ad.png'
#'After combining all outputs with an illustrator software, we get
#'Figure 1 (including `r fname_fig1` it if exists).
#+ include_fig1_illu
if (file.exists(fname_fig1)) {
  knitr::include_graphics(fname_fig1, dpi = 200)
}

#'## Short caption
#'
#'A) distribution of Error over unstacked Ridge learners. One can see that
#' Stacking  reduces the error substantially.
#'B) Breakdown of different stacking 
#' models. Note that models that combine power and connectivity perform best.
#'C) Variable importance for the full stacking model. The ridge models that
#' concatenate all frequencies were most important. However, single
#' frequency #'bands also contributed. The top features were  all related to power or
#' cross-power (covariance).
#'D) Partial dependency suggests monotonic but not linear relations between age
#' inputs and prediction. Note that some inputs have smaller value ranges,
#' which is consistent with data-driven regularization through the Ridge.
#'
#'## Conclusions so far.
#'
#'Let's revisit our questions
#'
#'1. What is the impact of stacking?
#'
#'Stacking greatly improves linear fits.
#'
#'2. Which type of MEG feature is most influential on the results?
#'
#'The Ridge model across frequency bands.
#'
#'3. Is the feature-type or the frequency band more impotant?
#'
#'There is substantial variability across frequencies and across feature types.
#'We see that beta and alpha band power and connectivity influenced the
#'performance.
#'
#'4. Does the stacker learn additive or interactive effects between features?
#'
#'The analysis of variable importance suggests some synergystic effects,
#'as features with low marginal performance make it among the top features.
#'
#'## Todos
#'
#'### global todo.
#'
#' - [ ] make plot of age prediction distributions for Ridge inputs
#' - [ ] think if we need more panels, e.g. to  show additional importance or to
#'  show a difference plot for models.
#' - [X] double check the extreme errors in sensor space models
#'    - It turns out scoring was based on input data, not age predictions.
#' - [ ] make really consistent x/y labels, e.g. input, output age  
#' - [ ] make sure, same units scatter plots have square layout.
#'
#'### Fig 1a stack models
#'
#'- [ ] overthink color code. Currently sky blue is double occupied.
#'
#'### Fig 1b stack models
#'
#'- [x] show essential, non-redundant models only, move rest to supplement
#'<!-- - [] add quantiles or mean info -->
#'- [ ] think of way to implement inference, visually, if needed supplement
#'- [x] simplify names
#'
#'### Fig 1c importance
#'
#'- [ ] choose which importances to show in main text
#'- [ ] choose model, or do inset or cut y axis
#'- [ ] produce elements for supplement
#'
#'### Figure 1d
#'- [x] add pdp plots
#'

#' ## Session info

#+ session_info
print(sessionInfo())
