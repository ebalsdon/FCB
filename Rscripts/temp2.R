# Run regression with no intercept
model3 <- lm(FTES ~ 0 + Frosh + Soph + Junior + Senior, data = data)

# View summary
summary(model3)

# Generate predicted FTES
data$FTES_pred3 <- predict(model3)

library(ggplot2)

ggplot(data, aes(x = Total)) +
  geom_point(aes(y = FTES), color = "blue", size = 2) +
  geom_point(aes(y = FTES_pred3), color = "purple", shape = 17, size = 2) +
  labs(title = "Observed vs. Predicted FTES (Model 3: by Student Level)",
       x = "Total",
       y = "FTES") +
  theme_minimal()
