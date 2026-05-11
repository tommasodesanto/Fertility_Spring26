# Plot annual birth probability by age and location type
library(ggplot2)
library(dplyr)

# Read the ASFR data
asfr <- read.csv("calibration_targets_output/asfr_by_location.csv")

# Convert to probability (divide by 1000) and clean up
asfr <- asfr %>%
  mutate(
    birth_prob = birth_rate / 1000 * 100,  # Convert to percentage
    location = case_when(
      location_type == "1_Peripheral" ~ "Peripheral",
      location_type == "2_Secondary" ~ "Secondary",
      location_type == "3_Superstar" ~ "Superstar"
    ),
    location = factor(location, levels = c("Peripheral", "Secondary", "Superstar")),
    # Create numeric midpoint for plotting
    age_mid = case_when(
      age_group_5yr == "15-19" ~ 17,
      age_group_5yr == "20-24" ~ 22,
      age_group_5yr == "25-29" ~ 27,
      age_group_5yr == "30-34" ~ 32,
      age_group_5yr == "35-39" ~ 37,
      age_group_5yr == "40-44" ~ 42,
      age_group_5yr == "45-50" ~ 47
    )
  ) %>%
  filter(age_group_5yr != "15-19" & age_group_5yr != "45-50")  # Focus on prime childbearing ages

# Create the plot
p <- ggplot(asfr, aes(x = age_mid, y = birth_prob, color = location, shape = location)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3.5) +
  scale_color_manual(
    values = c("Peripheral" = "#2E86AB", "Secondary" = "#A23B72", "Superstar" = "#F18F01"),
    name = "Location"
  ) +
  scale_shape_manual(
    values = c("Peripheral" = 16, "Secondary" = 17, "Superstar" = 15),
    name = "Location"
  ) +
  scale_x_continuous(
    breaks = c(22, 27, 32, 37, 42),
    labels = c("20-24", "25-29", "30-34", "35-39", "40-44")
  ) +
  scale_y_continuous(
    limits = c(0, 13),
    breaks = seq(0, 12, 2),
    labels = function(x) paste0(x, "%")
  ) +
  labs(
    x = "Age",
    y = "Annual probability of birth",
    title = "Fertility Timing by Location Type"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = c(0.85, 0.85),
    legend.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(color = "gray40", size = 11)
  )

# Save the plot
ggsave("calibration_targets_output/fertility_by_age_location.pdf", p, width = 8, height = 5)
ggsave("calibration_targets_output/fertility_by_age_location.png", p, width = 8, height = 5, dpi = 300)

cat("Figures saved to calibration_targets_output/\n")
print(p)
