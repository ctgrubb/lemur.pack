library(dplyr)

set.seed(66374)

N <- 1000
n <- 100

ex1_df1 <- data.frame(
  Group = sample(LETTERS[1:3], size = N, replace = TRUE, prob = c(0.35, 0.5, 0.15))
)

centers <- ifelse(ex1_df1$Group == "A", 120000, ifelse(ex1_df1$Group == "B", 240000, 400000))

ex1_df1$HousePrice <- rnorm(N, 0, 40000) + centers
ex1_df1$Income <- rnorm(N, mean = 0.2 * ex1_df1$HousePrice + 30000, sd = 20000)

# plot(ex1_df1$Income, ex1_df1$HousePrice)


ex1_df1 <- ex1_df1 %>%
  mutate(
    IncomeCat = case_when(
      Income < 40000 ~ "<40000",
      Income < 65000 ~ "40000-64999",
      Income < 110000 ~ "65000-109999",
      TRUE ~ ">=110000"
    ),
    HousePriceCat = case_when(
      HousePrice < 80000 ~ "<80000",
      HousePrice < 140000 ~ "80000-139999",
      HousePrice < 220000 ~ "140000-219999",
      HousePrice < 360000 ~ "220000-359999",
      TRUE ~ ">=360000"
    )
  ) %>%
  mutate(
    IncomeCat = factor(IncomeCat, levels = c("<40000", "40000-64999", "65000-109999", ">=110000")),
    HousePriceCat = factor(HousePriceCat, levels = c("<80000", "80000-139999", "140000-219999", "220000-359999", ">=360000"))
  ) %>% select(
    Group, Income, IncomeCat, HousePrice, HousePriceCat
  )

# table(ex1_df2)

pop <- ex1_df1 %>%
  arrange(HousePrice) %>%
  select(HousePrice, HousePriceCat)

samp <- ex1_df1 %>%
  sample_n(n) %>%
  select(
    Income,
    IncomeCat,
    HousePriceCat
  )

ex1samp <- samp
ex1pop <- pop

usethis::use_data(ex1samp, overwrite = TRUE)
usethis::use_data(ex1pop, overwrite = TRUE)
