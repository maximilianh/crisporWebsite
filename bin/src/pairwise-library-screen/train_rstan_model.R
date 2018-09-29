#  
# Usage: Rscript train_rstan_model.R input_data.txt guide_length param_values.txt
# 
# The input data file should be prepared using encode_full_matrix.py
# guide_length specifies the guide length that will be used for the model
#


library(rstan)
library(plyr)
library(reshape2)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


 # Parse arguments

args <- commandArgs(TRUE)

full.df <- read.table(args[1])
guide.length.filter <- as.integer(args[2])
out.file <- args[3]

num.cores <- parallel::detectCores()

# Run model

full.df <- full.df[complete.cases(full.df),]

names(full.df)[1:7] <- c('Indel.Rate', 'Target.Ratio', 'Guide.Group', 'Guide.Length', 'Guide.Seq', 'Target.Seq', 'Num.MM')

meta.df <- full.df[,1:7]

data.df <- full.df[meta.df$Guide.Length == guide.length.filter, 8:ncol(full.df)]
meta.df <- meta.df[meta.df$Guide.Length == guide.length.filter,]

meta.df$Guide.Group <- meta.df$Guide.Group + 1

X <- as.matrix(data.df)

N <- nrow(meta.df)

L <- guide.length.filter

model.data <- list(N = N, L = L, 
				   num_guides = max(meta.df$Guide.Group),
				   matrix_cols = ncol(X),
				   guide_index = meta.df$Guide.Group,
				   indel_rate = meta.df$Indel.Rate,
				   mm_matrix = X,
				   beta_hyperparameter = 1,
				   sigma = 0.01)

stan.fit <- stan(file = "exp_model.stan", data = model.data, iter = 2500, chains = num.cores,
  thin = 2, warmup = 500,  pars = c("mm_penalty", "guide_effect", "pred_indel", 'mult_effect'))

save.image(paste0('stan_temp_', guide.length.filter,'.Rdata'))

rev.pos <- guide.length.filter:1

stan.summ <- summary(stan.fit)

mm.index <- c('dT:rC', 'dT:rG', 'dT:rU', 'dG:rA', 'dG:rG', 'dG:rU', 'dC:rA', 'dC:rC', 'dC:rU', 'dA:rA', 'dA:rC', 'dA:rG')


mult.ci.df <- as.data.frame(stan.summ$summary[grepl('mult', rownames(stan.summ$summary)), c('mean', 'se_mean', '2.5%', '50%', '97.5%')])
mult.ci.df$Mismatch <- mm.index
mult.ci.df$Position <- 0:(nrow(mult.ci.df) - 1) %/% 12 + 1
mult.ci.df <- subset(mult.ci.df, Position <= guide.length.filter)
mult.ci.df$Position <- rev.pos[mult.ci.df$Position]
names(mult.ci.df)[3:5] <- c('Low', 'Median', 'High')

write.table(mult.ci.df, out.file, quote = FALSE, sep = '\t', row.names = FALSE)


