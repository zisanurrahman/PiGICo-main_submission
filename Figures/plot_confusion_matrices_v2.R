#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(ggplot2); library(tidyr); library(optparse)
})

opt <- list(
  make_option(c("-p","--preds"), type="character", help="predictions.csv"),
  make_option(c("-t","--task"),  type="character", default="HEALTH", help="Task name"),
  make_option(c("-o","--out"),   type="character", default="confusion_mats", help="Output prefix"),
  make_option(c("--threshold"),  type="double", default=NA_real_,
              help="Decision threshold (default = Youden J per model over CV)")
)
opt <- parse_args(OptionParser(option_list=opt))
stopifnot(!is.null(opt$preds) && file.exists(opt$preds))

pred <- read_csv(opt$preds, show_col_types = FALSE) %>% filter(Task == opt$task)
models <- c("LR","RF","XGB"); pred <- pred %>% filter(Model %in% models)

# Pick threshold
choose_thresh <- function(df) {
  if (!is.na(opt$threshold)) return(opt$threshold)
  # Youden J on pooled folds per model
  grid <- seq(0,1,by=0.001)
  sapply(split(df, df$Model), function(d) {
    y <- d$y_true; p <- d$y_prob
    J <- sapply(grid, function(th) {
      pr <- as.integer(p >= th)
      TP <- sum(pr==1 & y==1); FP <- sum(pr==1 & y==0)
      TN <- sum(pr==0 & y==0); FN <- sum(pr==0 & y==1)
      TPR <- ifelse((TP+FN)>0, TP/(TP+FN), 0)
      FPR <- ifelse((FP+TN)>0, FP/(FP+TN), 1)
      TPR - FPR
    })
    grid[which.max(J)]
  })
}
thresh_by_model <- choose_thresh(pred)

# Confusion tables per model (pooled over folds)
cm_list <- lapply(models, function(m) {
  d <- pred %>% filter(Model==m)
  th <- ifelse(is.na(opt$threshold), thresh_by_model[[m]], opt$threshold)
  d <- d %>% mutate(y_hat = as.integer(y_prob >= th))
  TP <- sum(d$y_hat==1 & d$y_true==1)
  FP <- sum(d$y_hat==1 & d$y_true==0)
  TN <- sum(d$y_hat==0 & d$y_true==0)
  FN <- sum(d$y_hat==0 & d$y_true==1)
  acc <- (TP+TN)/nrow(d)
  tibble(Model=m, TN=TN, FP=FP, FN=FN, TP=TP, ACC=acc, TH=th)
})
cm <- bind_rows(cm_list)

# Heatmap style faceted confusion matrices
plot_df <- cm %>%
  pivot_longer(cols=c("TN","FP","FN","TP"), names_to="Cell", values_to="N") %>%
  mutate(Row = ifelse(Cell %in% c("TN","FN"), "True 0", "True 1"),
         Col = ifelse(Cell %in% c("TN","FP"), "Pred 0", "Pred 1"))

p <- ggplot(plot_df, aes(x=Col, y=Row, fill=N)) +
  geom_tile(color="white", linewidth=0.6) +
  geom_text(aes(label=N), size=5) +
  scale_fill_gradient(low="#eef2ff", high="#1B9E77") +
  facet_wrap(~Model, nrow=1) +
  coord_equal() +
  labs(title=sprintf("Confusion matrices — %s (pooled CV)", opt$task),
       x=NULL, y=NULL, fill="Count") +
  theme_minimal(base_size=12) +
  theme(axis.text = element_text(color = 'black', size= 14))+
  theme(panel.grid=element_blank(), strip.text=element_text(face="bold"))

ggsave(paste0(opt$out,"_",opt$task,".png"), p, width=8, height=3, dpi=600,bg='white')
ggsave(paste0(opt$out,"_",opt$task,".pdf"), p, width=10, height=3.2)

