# ==========================================
# پروژه تحلیل چندمتغیره - ASD Dataset
# نسخه اصلاح‌شده (بدون متغیر result که باعث سینگولار شدن کوواریانس می‌شد)
# ==========================================

library(ggplot2)
library(corrplot)
library(MASS)
library(psych)
library(GGally)
library(mvnormtest)
library(biotools)

# ==========================================
# لود و پاکسازی داده
# ==========================================

asd <- read.csv("https://raw.githubusercontent.com/shaheennamboori/CSV_dataset_for_autism_diagnostics/master/Autism-Adult_Data.csv")

asd$age[asd$age == "?"] <- NA
asd$age <- as.numeric(asd$age)
asd$Class.ASD <- as.factor(asd$Class.ASD)
asd_clean <- asd[!is.na(asd$age), ]

# متغیرهای کمی (result حذف شد چون result = sum(A1..A10) و باعث هم‌خطی کامل می‌شد)
num_vars <- asd_clean[, c("A1_Score","A2_Score","A3_Score","A4_Score","A5_Score",
                          "A6_Score","A7_Score","A8_Score","A9_Score","A10_Score",
                          "age")]

n <- nrow(num_vars)
p <- ncol(num_vars)

# ==========================================
# الف) تحلیل توصیفی چندمتغیره
# ==========================================

cat("\n=== بردار میانگین ===\n")
print(round(colMeans(num_vars), 3))

cat("\n=== ماتریس کوواریانس ===\n")
cov_matrix <- cov(num_vars)
print(round(cov_matrix, 3))

cat("\n=== ماتریس همبستگی ===\n")
cor_matrix <- cor(num_vars)
print(round(cor_matrix, 3))

corrplot(cor_matrix, method = "color", type = "upper",
         title = "Correlation Matrix - ASD Dataset", mar = c(0,0,2,0))

pairs(num_vars[,1:6],
      col = ifelse(asd_clean$Class.ASD == "YES", "red", "blue"),
      main = "Scatterplot Matrix - ASD Scores")

cat("\n=== داده‌های پرت (فاصله ماهالانوبیس) ===\n")
mah_dist <- mahalanobis(num_vars, colMeans(num_vars), cov(num_vars))
cutoff <- qchisq(0.975, df = p)
outliers <- which(mah_dist > cutoff)
cat("تعداد داده‌های پرت:", length(outliers), "\n")
cat("ردیف‌های پرت:", outliers, "\n")

plot(mah_dist, type = "p", pch = 19, cex = 0.5,
     main = "Mahalanobis Distance - Outlier Detection",
     ylab = "Mahalanobis Distance", xlab = "Observation")
abline(h = cutoff, col = "red", lty = 2)
legend("topright", legend = "Cutoff (97.5%)", col = "red", lty = 2)

# ==========================================
# ب) نرمال چندمتغیره (Mardia دستی)
# ==========================================

S <- cov(num_vars)
X_c <- scale(num_vars, center = TRUE, scale = FALSE)
A_mat <- X_c %*% solve(S) %*% t(X_c)

b1p <- sum(A_mat^3) / n^2
b2p <- sum(diag(A_mat^2)) / n
k_skew <- n * b1p / 6
df_skew <- p * (p+1) * (p+2) / 6
p_skew <- pchisq(k_skew, df = df_skew, lower.tail = FALSE)
z_kurt <- (b2p - p*(p+2)) / sqrt(8*p*(p+2)/n)
p_kurt <- 2 * pnorm(abs(z_kurt), lower.tail = FALSE)

cat("\n=== آزمون نرمال چندمتغیره (Mardia) ===\n")
cat("Mardia Skewness:", round(b1p,3), "| Chi-sq:", round(k_skew,3), "| p-value:", round(p_skew,5), "\n")
cat("Mardia Kurtosis:", round(b2p,3), "| z:", round(z_kurt,3), "| p-value:", round(p_kurt,5), "\n")
if (p_skew > 0.05 & p_kurt > 0.05) {
  cat("نتیجه: داده‌ها نرمال چندمتغیره هستند\n")
} else {
  cat("نتیجه: داده‌ها نرمال چندمتغیره نیستند\n")
}

qqplot(qchisq(ppoints(n), df = p), mah_dist,
       main = "Chi-Square Q-Q Plot (Multivariate Normality)",
       xlab = "Theoretical Quantiles", ylab = "Mahalanobis Distance")
abline(0, 1, col = "red")

# ==========================================
# ج) آزمون فرض چندمتغیره
# ==========================================

cat("\n=== آزمون Hotelling T2 ===\n")
group1 <- num_vars[asd_clean$Class.ASD == "YES", ]
group2 <- num_vars[asd_clean$Class.ASD == "NO", ]
n1 <- nrow(group1); n2 <- nrow(group2)
mean1 <- colMeans(group1); mean2 <- colMeans(group2)
S1 <- cov(group1); S2 <- cov(group2)
Sp <- ((n1-1)*S1 + (n2-1)*S2) / (n1+n2-2)
diff <- mean1 - mean2
T2 <- as.numeric((n1*n2/(n1+n2)) * t(diff) %*% solve(Sp) %*% diff)
F_stat <- T2 * (n1+n2-p-1) / ((n1+n2-2)*p)
df1 <- p; df2 <- n1+n2-p-1
p_value <- pf(F_stat, df1, df2, lower.tail = FALSE)

cat("T2 =", round(T2, 3), "\n")
cat("F =", round(F_stat, 3), "\n")
cat("df1 =", df1, ", df2 =", df2, "\n")
cat("p-value =", format.pval(p_value, digits = 4), "\n")
if (p_value < 0.05) cat("نتیجه: تفاوت معناداری بین دو گروه وجود دارد\n")

cat("\n=== MANOVA ===\n")
score_matrix <- as.matrix(asd_clean[, c("A1_Score","A2_Score","A3_Score",
                                        "A4_Score","A5_Score","A6_Score",
                                        "A7_Score","A8_Score","A9_Score",
                                        "A10_Score","age")])
manova_result <- manova(score_matrix ~ Class.ASD, data = asd_clean)
print(summary(manova_result, test = "Wilks"))
print(summary(manova_result, test = "Pillai"))

# ==========================================
# د) تحلیل مولفه‌های اصلی PCA
# ==========================================

cat("\n=== PCA ===\n")
pca_result <- prcomp(num_vars, scale. = TRUE)
print(summary(pca_result))

eigenvalues <- pca_result$sdev^2
plot(eigenvalues, type = "b", pch = 19,
     main = "Scree Plot - ASD Dataset",
     xlab = "Principal Component", ylab = "Eigenvalue")
abline(h = 1, col = "red", lty = 2)
legend("topright", legend = "Eigenvalue = 1", col = "red", lty = 2)

biplot(pca_result, main = "PCA Biplot - ASD Dataset", cex = 0.6)

cat("\n=== بارگذاری مولفه‌ها (4 مولفه اول) ===\n")
print(round(pca_result$rotation[,1:4], 3))

# ==========================================
# هـ) تحلیل عاملی
# ==========================================

cat("\n=== تحلیل عاملی ===\n")
fa_result <- fa(num_vars, nfactors = 3, rotate = "varimax", fm = "ml")
print(fa_result$loadings, cutoff = 0.3)
fa.diagram(fa_result, main = "Factor Analysis - ASD Dataset")

# ==========================================
# ز) تحلیل ممیزی
# ==========================================

cat("\n=== تحلیل ممیزی (LDA) ===\n")
lda_result <- lda(Class.ASD ~ A1_Score+A2_Score+A3_Score+A4_Score+A5_Score+
                    A6_Score+A7_Score+A8_Score+A9_Score+A10_Score+age,
                  data = asd_clean)
print(lda_result)

lda_pred <- predict(lda_result)
conf_matrix <- table(Actual = asd_clean$Class.ASD, Predicted = lda_pred$class)
cat("\n=== ماتریس طبقه‌بندی ===\n")
print(conf_matrix)

accuracy <- sum(diag(conf_matrix)) / sum(conf_matrix) * 100
cat("\nدرصد صحت طبقه‌بندی:", round(accuracy, 2), "%\n")

lda_df <- data.frame(LD1 = lda_pred$x[,1], Class = asd_clean$Class.ASD)
ggplot(lda_df, aes(x = LD1, fill = Class)) +
  geom_histogram(bins = 30, alpha = 0.6, position = "identity") +
  scale_fill_manual(values = c("NO" = "blue", "YES" = "red")) +
  labs(title = "LDA - توزیع تابع ممیزی", x = "LD1", y = "فراوانی") +
  theme_minimal()

cat("\n✓ تحلیل کامل شد!\n")
