# Load necessary libraries
library(ggplot2)
library(dplyr)

# Load the data (assuming the file is named 'carnivoreteeth.csv' and is in your working directory)
carnivore_data <- read.csv("carnivoreteeth.csv", header= TRUE)
head(carnivore_data)
colnames(carnivore_data)[1] <- "Species"
data0<-subset(carnivore_data, Lat>0)

# Create the scatter plot
fig1<-ggplot(data0, aes(x = Lat, y = PM4, color = Species)) +
  geom_point(size = 3, alpha = 0.7) +  # Scatter points
  geom_smooth(aes(group = Species), method = "lm", se = FALSE, linetype = "dashed") +  # Trend lines for each species
  geom_smooth(method = "lm", se = FALSE, color = "black", size = 1) +  # Common trend line
  labs(title = "Relationship between Latitude and PM4",
       x = "Latitude",
       y = "PM4",
       color = "Species") +
  theme_minimal() +
  theme(legend.position = "right")
tiff(filename = "initial_model.tiff", width = 30, height = 10, units = "in", res = 300)
fig1
dev.off()


data1<- subset(data0, sex=male)
# Create the scatter plot
fig2<-ggplot(data1, aes(x = PM4, y = CsupL, color = Family)) +
  geom_point(size = 3, alpha = 0.7) +  # Scatter points
  geom_smooth(aes(group = Family), method = "lm", se = FALSE, linetype = "dashed") +  # Trend lines for each species
  geom_smooth(method = "lm", se = FALSE, color = "black", size = 1) +  # Common trend line
  labs(title = "Relationship between CsupL and PM4",
       x = "PM4",
       y = "CsupL",
       color = "Family") +
  theme_minimal() +
  theme(legend.position = "right")
tiff(filename = "initial_model_CsupL_PM4.tiff", width = 30, height = 10, units = "in", res = 300)
fig2
dev.off()
