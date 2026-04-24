######################################################################################
# Predict RMR, MMR, and aerobic scope using temperature and fish size
# Extrapolate to 1) warmer temperatures, and 2) larger fish
# Demer and Zwolinski 2014 show that mean landed sardine lengths off San Pedro ~ 14cm
# Using length-weight relationship in Kuriyama et al. 2024 assessment:
# tl.cm <- (3.574 + SL.mm * 1.149) / 10 = 16.4434 cm TL
# and wgt(kg) = 4.446313e-06 * (tl.cm ^ 3.197) = 0.034318 kg = 34.3 g
# Contact bmuhling@ucsc.edu
######################################################################################

library(readxl)
library(mgcv)
library(ggplot2)
library(viridis)
library(gratia)

# Load and combine the data (data are not public)
t135 <- read_xlsx("./caClimate/data/lonthair/1_OutputFiles/Sardine_Energetics_060522.xlsx", sheet = "13.5")
t165 <- read_xlsx("./caClimate/data/lonthair/1_OutputFiles/Sardine_Energetics_060522.xlsx", sheet = "16.5")
t195 <- read_xlsx("./caClimate/data/lonthair/1_OutputFiles/Sardine_Energetics_060522.xlsx", sheet = "19.5")
# Use shorter colnames 
colnames(t135) <- colnames(t165)
t <- rbind(t135, t165, t195)
# Better colnames for fields we need
t$rmr <- t$`Individual Absolute RMR`
t$wgtRMR <- t$`Ind Weight (RMR)`
t$aas <- t$AAS
t$mmr <- t$`Individual Absolute MMR - 1 minute`

# Predict RMR using temp and fish size (mg O2 / hr per fish) using a GAM
gamRMR <- gam(rmr ~ s(wgtRMR) + s(Temp, k = 3), data = t, family = Gamma(link = "log"))
summary(gamRMR) # dev expl 95.6%
plot(gamRMR, scale = 0, pages = 1) 
gratia::draw(gamRMR) & theme_bw() 
# Save plot
ggsave(filename = "./caClimate/plots/gamRMR_partials.tiff", dpi = 600, compression = "lzw",
       width = 6, height = 4)

# Predict MMR using temp and fish size (mg O2 / hr per fish)
gamMMR <- gam(mmr ~ s(wgtRMR) + s(Temp, k = 3), data = t, family = Gamma(link = "log"))
summary(gamMMR) # 91.7%
plot(gamMMR, scale = 0, pages = 1) 
gratia::draw(gamMMR) & theme_bw() 
# Save plot
ggsave(filename = "./caClimate/plots/gamMMR_partials.tiff", dpi = 600, compression = "lzw",
       width = 6, height = 4)

# Save models
save(gamRMR, file = "./caClimate/models/gamRMR.rda")
save(gamMMR, file = "./caClimate/models/gamMMR.rda")

# Show what happens if we extrapolate temperature
toScore <- data.frame(expand.grid("wgtRMR" = seq(20, 120, by = 20), "Temp" = seq(12, 28, by = 2)))
toScore$predRMR <- predict(gamRMR, toScore, type = "response")
toScore$predMMR <- predict(gamMMR, toScore, type = "response")
toScore$AAS <- toScore$predMMR - toScore$predRMR
# AAS as a % of RMR
toScore$AASpercRMR <- toScore$AAS / toScore$predRMR
# AAS as a % of MMR
toScore$AASpercMMR <- toScore$AAS / toScore$predMMR
# Show RMR as a % of max scope
toScore$percRMR <- toScore$predRMR / toScore$predMMR
  
# Predicted RMR
ggplot(toScore) + geom_tile(aes(x = wgtRMR, y = Temp, fill = predRMR)) + 
  scale_fill_viridis("Predicted \nRMR") + theme_bw() + xlab("Fish Weight (g)") + ylab("Water Temperature (°C)") + 
  ggtitle("Predictions of routine metabolic rate") + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  theme(plot.title = element_text(size = 14))
# Save plot
ggsave(filename = "./caClimate/plots/pred_rmr.tiff", dpi = 600, compression = "lzw",
       width = 4, height = 3)

# Predicted MMR
ggplot(toScore) + geom_tile(aes(x = wgtRMR, y = Temp, fill = predMMR)) + 
  scale_fill_viridis("Predicted \nMMR") + theme_bw() + xlab("Fish Weight (g)") + ylab("Water Temperature (°C)") +
  ggtitle("Predictions of maximum metabolic rate") + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  theme(plot.title = element_text(size = 14))
# Save plot
ggsave(filename = "./caClimate/plots/pred_mmr.tiff", dpi = 600, compression = "lzw",
       width = 4, height = 3)

# Predicted RMR as a % of MMR (eg 0.2 means they're using 20% of their total scope
ggplot(toScore) + geom_tile(aes(x = wgtRMR, y = Temp, fill = percRMR)) + 
  scale_fill_viridis("", breaks = c(0.08, 0.23), labels = c("Higher \nMetabolic \nScope", "Lower \nMetabolic \nScope")) + 
  theme_bw() + xlab("Fish Weight (g)") + ylab("Water Temperature (°C)") + 
  ggtitle("Routine metabolic rate as % of \nmaximum metabolic rate") + 
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  theme(plot.title = element_text(size = 14))
# Save plot
ggsave(filename = "./caClimate/plots/pred_met_scope.tiff", dpi = 600, compression = "lzw",
       width = 4, height = 3)



