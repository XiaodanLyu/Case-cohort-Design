## developed by @XiaodanLyu
sample <- read.csv("data/realdata.csv", header = T)
sample$delta <- factor(sample$delta, levels = c(1, 0), labels = c("relapse", "censored"))
names(sample)[2] <- "Observation"
require(ggplot2)
ggplot(data = sample[1:100, ], aes(x = age, y = Time, pch = Observation, col = Observation)) +
  geom_point(cex = 3) + xlab("Age (year)") + ylab("Relapse Time (day)\n") +
  theme(text = element_text(size=20),axis.title.x=element_text(vjust=-0.25)) 
ggsave("Sample.png", width = 8.5, height = 5.5)

source("functions/cc-estimation.r")
source("functions/cox-estimation.r")
X <- read.csv("data/realdata.csv", header = T)
names(X)[1] <- "time"
cc.est <- cc.estimation(X, c(0, 0, 0))
Y <- read.csv("data/coxsample.csv", header = T)
cox.est <- cox.estimation(Y, c(0, 0, 0))
cc.est2 <- cc.estimation(X[, -6], c(0, 0))
cox.est2 <- cox.estimation(Y[, -5], c(0, 0))
