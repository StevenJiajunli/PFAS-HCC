
TCGA_tri_OS <- cbind(TCGA_ML,predictions)

dat <- TCGA_tri_OS

# Sort the data frame by the 'cox_score' column in ascending order
dat <- dat[order(dat$`1`), ]

dat$Samplerisk <- c(1:nrow(dat))
dat$riskScore <- dat$`1`
cutpoint <- round(median(dat$riskScore),5)
dat$risk_group[dat$riskScore >= cutpoint] <- "High PFASHRSig"
dat$risk_group[dat$riskScore < cutpoint] <- "Low PFASHRSig"

p1 <- ggplot(dat, aes(Samplerisk, riskScore)) +
  # Plot the points, filling based on risk group
  geom_point(aes(fill = risk_group), pch = 21, color = "white", size = 3, stroke = 0.1) +
  
  # Set fill colors for risk groups
  scale_fill_manual(values = c("#eea9b6","#44a5cb")) +
  
  # Label axes
  labs(x = "Patient ID (increasing Risk Score)", y = "Risk Score") +
  
  # Annotate cutpoint value
  annotate("text", x = sum(dat$risk_group == "Low PFASHRSig"), y = cutpoint + 0.2, size = 5.5,
           label = paste("cutpoint=", cutpoint)) +
  
  # Add horizontal line at cutpoint
  geom_hline(yintercept = cutpoint, colour = "black", linetype = "dashed", size = 0.8) +
  
  # Add vertical line at cutpoint
  geom_vline(xintercept = sum(dat$risk_group == "Low PFASHRSig"), colour = "black", linetype = "dashed", size = 0.8) +
  
  # Customize y-axis
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.01)),
                     breaks = seq(-2, max(dat$riskScore) + 0.5, 0.5)) +
  
  # Customize x-axis
  scale_x_continuous(expand = expansion(mult = c(0.03, 0.03)),
                     limits = c(0, max(dat$Samplerisk) + 1)) +
  
  # Use white background and adjust theme
  theme_bw() +
  theme(
    axis.title = element_text(size = 13, color = 'black'),
    axis.text = element_text(size = 10, color = 'black'),
    legend.title = element_blank(),
    legend.text = element_text(size = 12, color = 'black'),
    legend.background = element_blank(),
    legend.position = c(0.2, 0.85)
  ) +
  
  # Adjust legend settings
  guides(fill = guide_legend(ncol = 1, override.aes = list(size = 3)))

# Print the plot
p1

# Create a new column 'event' based on the 'fustat' column
dat$event <- ifelse(dat$OS == "1", "death", "alive")

# Create the plot p2
p2 <- ggplot(dat, aes(x = Samplerisk, y = OS.time)) +
  
  # Add points with color and shape based on 'event' (death or alive)
  geom_point(aes(fill = event, shape = event), color = "white", size = 3.5, stroke = 1) +
  
  # Set custom fill colors for 'event' categories
  scale_fill_manual(values = c("#44a5cb", "#eea9b6")) +
  
  # Set custom shapes for 'event' categories
  scale_shape_manual(values = c(22, 21)) +
  
  # Add labels for axes
  labs(x = "Patient ID (increasing risk score)", y = "Survival time (Month)") +
  
  # Add a vertical line at the point where risk group changes (low risk)
  geom_vline(xintercept = sum(dat$risk_group == "Low PFASHRSig"), colour = "black", linetype = "dashed", size = 0.8) +
  
  # Customize x-axis
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01)),
                     limits = c(0, max(dat$Samplerisk) + 1)) +
  
  # Use a white background theme
  theme_bw() +
  
  # Customize plot theme
  theme(
    axis.title = element_text(size = 13, color = 'black'),
    axis.text = element_text(size = 10, color = 'black'),
    legend.title = element_blank(),
    legend.text = element_text(size = 12, color = 'black'),
    legend.background = element_blank(),
    legend.position = c(0.15, 0.85)
  ) +
  
  # Customize legend
  guides(fill = guide_legend(ncol = 1, override.aes = list(size = 5)))

# Print the plot
p2

library(ggh4x)
library(tidyr)
library(dplyr)
library(ggplot2)

final_sig <- c("ESR1", "APOA1", "IGF1", "PPARGC1A", "SERPINE1", 
               "PON1", "HMOX1", "APCS", "ACADS", "SLC10A1", "SLC2A2", "ACAT1", "C1S", "LCAT")

dat <- dat[,c("OS.time","OS",final_sig,"Samplerisk","riskScore","risk_group","event")]
colnames(dat)

# Prepare the heat data
# Standardize the data (z-scores)
heat <- dat[, 3:16]%>%
  scale() %>%                        # Step 1: Standardize the data
  as.data.frame()                     # Step 2: Convert the scaled matrix to a data frame

# Add row names as a column called 'Sample'
heat <- heat %>%
  rownames_to_column("Sample")   # Step 3: Extract row names as 'Sample'

# Convert 'Sample' to a factor
heat <- heat %>%
  mutate(Sample = factor(Sample, levels = Sample)) # Step 4: Convert 'Sample' to a factor

# Reshape the data into long format (key = "gene", value = "expression")
heat <- heat %>%
  gather(key = "gene", value = "expression", -Sample)  # Step 5: Reshape the data

dat$Sample <- rownames(dat)

# Perform an inner join with 'dat' to include the 'risk_group' column
heat <- heat %>%
  inner_join(., dat[, c("Sample", "risk_group")], by = c("Sample")) # Step 6: Merge with 'dat'

unique(dat$risk_group)
# Convert 'risk_group' to a factor and set levels
heat <- heat %>%
  mutate(risk_group = factor(risk_group, levels = unique(risk_group)))  # Step 7: Convert 'risk_group' to a factor

# Create the heatmap plot
p_heat <- ggplot(heat, aes(Sample, gene, fill = expression)) +
  geom_tile() +  # Create the heatmap tiles
  scale_fill_gradientn(colors = c('#4E72B8', '#ADD8E6', '#FFFFFF', '#FFA500', '#FF0000')) +  # Custom color scale
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.text.x.bottom = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") +
  xlab("") + ylab("") +
  theme(
    axis.title = element_text(size = 10, color = 'black'),
    axis.text = element_text(size = 9, color = 'black')
  ) +
  facet_nested(. ~ risk_group, drop = TRUE, scale = "free", space = "free", 
               strip = strip_nested(
                 background_x = elem_list_rect(fill = c("#44a5cb", "#eea9b6")),  # Different background colors
                 text_x = element_text(size = 10),
                 by_layer_x = FALSE
               ))

# Display the heatmap
p_heat

library(patchwork)
p1 / p2 / p_heat
