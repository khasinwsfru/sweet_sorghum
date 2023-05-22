library(car)
library(truncnorm)
library(mgcv)
library(tidyverse)
library(ggplot2)

myenzgam3 <- gam(lesion.mm~
                   factor(Assay)+
                   factor(Timepoint)+
                   s(umoltgcl_ugprotein_min)+
                   s(abs_abts_consumed)+
                   abs__H2O2used_ugprotein_min__+
                   s(nmolcinna_ug_min),
                 data = noday0, #i filtered out day0s otherwise predict() complains
                 family = gaussian(link="log")
                   )

plot(enzgam3)
gam.check(enzgam3)
anova(myenzgam3)


#make a whole vector of fake data - only use one row of fake data at a time though. pick a value for everything else.

myenztest  <- data.frame(Assay = c(rep("Assay1", 8), rep("Assay2", 9), rep("Assay1", 9), rep("Assay2", 8)),
                    Timepoint = c(rep("Expansion", 17), rep("Initation", 17)),
                    umoltgcl_ugprotein_min = rtruncnorm(enz3model$umoltgcl_ugprotein_min, a = 1, b = 10, mean = mean(enz3model$umoltgcl_ugprotein_min), sd = 3.7),
                    abs_abts_consumed = rtruncnorm(enz3model$abs_abts_consumed, a = 0.025, b = 0.25, mean(enz3model$abs_abts_consumed), sd = 0.2),
                    nmolcinna_ug_min = rtruncnorm(enz3model$nmolcinna_ug_min, a = 0.0025, b = 0.025, mean = mean(enz3model$nmolcinna_ug_min), sd = 0.006))


#this worked best but I tested others

paltest_prx_hi_abts_lo <- data.frame(Assay = c(rep("Assay1", 8), rep("Assay2", 9), rep("Assay1", 9), rep("Assay2", 8)),
                    Timepoint = c(rep("Expansion", 17), rep("Initation", 17)),
                    umoltgcl_ugprotein_min = 6,
                    abs_abts_consumed = 0.15,
                    nmolcinna_ug_min = myenztest$nmolcinna_ug_min)

myenzfits2 <- predict(myenzgam3, newdata=paltest_prx_hi_abts_lo, type='response', se=TRUE)

enzpredicts2 <- data.frame(paltest_prx_hi_abts_lo, myenzfits2) %>% 
  mutate(lower = fit - 1.96*se.fit,
         upper = fit + 1.96*se.fit)

  palplot <- ggplot(aes(x=nmolcinna_ug_min,y=fit), data=enzpredicts2) +
  geom_ribbon(aes(ymin = lower, ymax=upper), fill='gray90') +
  geom_line(color='#00aaff') +theme_bw()+
  facet_wrap(c("Assay", "Timepoint"))
