library(brms)
library(tidyverse)


# # Fitting Curves
# print("Preparing Data")

# getwd()

#Read Data
df <- read_csv("./data/Francisca//GrowthRate20C/Pairwise_M9_Experiment.csv", col_type = cols()) %>%
    select(!contains("Blank")) %>%
    pivot_longer(-Reading, "ID", "value") %>%
    separate(ID,into = c("OTU","Treatment","Replicate"), sep = "_")

#Get ID
Treat <- unique(df$Treatment)
ID <- unique(df$OTU)
Rep <- unique(df$Replicate)

fit_IDs <- expand_grid(Treat,ID, Rep)

#set priors
p <- prior(normal(0.1,0.1),nlpar = "b1", lb = 0)  +
     prior(normal(0.2,0.1),nlpar = "b2", lb = 0) +
     prior(normal(0.2,0.1), nlpar = "b3", lb = 0)

print("First Fit")
##
# Fit
##
fits <- list()

#first model
df_fit_first <- df %>%
    filter(OTU == fit_IDs$ID[1], Treatment == fit_IDs$Treat[1], Replicate == fit_IDs$Rep[1])

fits[[1]] <- brm(bf(log(value) ~ log((b1 * b2) / (b1 + (b2 - b1) * exp(-b3 * Reading))), b1 + b2 + b3 ~ 1, nl = TRUE),
           data = df_fit_first, prior = p, refresh = 0)


print("Second fit")
#Update for others
for(i in 2:nrow(fit_IDs)){
print(paste("fitting: ", i, " Of ", nrow(fit_IDs)))
    
   df_fit_new <- df %>%
       filter(OTU == fit_IDs$ID[i], Treatment == fit_IDs$Treat[i], Replicate == fit_IDs$Rep[i])
    
   fits[[i]] <- update(fits[[1]],newdata = df_fit_new, refresh = 0)
    
    
}

# print("Saving Fits")
saveRDS(fits, "./data/Fits/bayes_curves.RDS")

# fits <- readRDS("./data/Fits/bayes_curves.RDS")

#plots

for(i in 1:nrow(fit_IDs)){
    name <- paste(fit_IDs[i,1],fit_IDs[i,2],fit_IDs[i,3],sep = "_")
    name <- str_replace(name,"/","_")

    p1 <- plot(fits[[i]])[[1]]
    ggsave(paste0("./data/Fits/Sample_figs/",name,".pdf"), p1)

    p2 <- plot(conditional_effects(fits[[i]]), points = T)[[1]]
    ggsave(paste0("./data/Fits/Effects_figs/",name,".pdf"), p2)

    p3 <- plot(pp_check(fits[[i]]), points = T)
    ggsave(paste0("./data/Fits/pp_checks/",name,".pdf"), p3)    
}



# #total_fit
# x <- read_csv("Data/Interactions_data.csv")

# fit <- brm(bf(value ~ Treatment + (Treatment | pairs)), data = x, family = gaussian(),
#             cores = 4,control = list(adapt_delta = 0.99))

# saveRDS(fit,"Data/bayes_fit.RDS")


# #pair_level_fits
# fit_pair <- df_mod %>%
#     nest(-pairs) %>%
#     mutate(fit = map(data, ~brm(bf(value ~ Treatment), data = ., family = gaussian(),cores = 4)) ) 

# saveRDS(fit_pair, "Data/bayes_pairs.RDS")