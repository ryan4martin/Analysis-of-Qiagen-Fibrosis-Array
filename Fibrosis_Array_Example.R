# Required libraries
library(tidyverse)
library(readxl)
library(here)
library(ggpubr)

# Read in Ct values for all treatments
# Each row consists of the Ct values for a single gene (name is in first column)
# Each column is labelled with the treatment and biological replicate number
dat <- read_xlsx(here("test_dat.xlsx"), col_names = T) %>% as_tibble()

# Convert all columns, except for first column with gene names, to numeric
dat[,2:ncol(dat)] <- apply(dat[,2:ncol(dat)], 2, as.numeric)

# Drop rows with NA in any of the columns
dat <- dat %>% drop_na()

# Select rows with all values less than 35 
dat <- dat %>% filter_if(is.numeric, .vars_predicate =  any_vars(abs(.) < 35))

# Pivot table to long format and split column with treatment time, 
# conditions and biological replicate into separate rows
# Remove control rows and all potential housekeeping genes from dataset
FibrosisArrayCt <- dat %>% pivot_longer(-Gene, names_to = "treatment", values_to = "Ct") %>%
  separate(treatment, into = c("Time", "Treatment", "n"), sep = "_") %>%
  separate(Treatment, into = c("Agonist", "siRNA"), sep = "/") %>%
  filter(!(Gene %in% c("RTC", "PPC", "Actb", "B2m", "Rplp1", 'Hprt1', 'Ldha')))

# Process dataset as before, this time select rows with housekeeping genes used for analysis
# Find mean Ct of housekeeping genes within each treatment group
SummaryHouseKeeping <- dat %>% pivot_longer(-Gene, names_to = "treatment", values_to = "Ct") %>%
  separate(treatment, into = c("Time", "Treatment", "n"), sep = "_") %>%
  separate(Treatment, into = c("Agonist", "siRNA"), sep = "/") %>%
  filter(Gene %in% c('Hprt1', 'Ldha'))  %>%
  group_by(Time, n, Agonist, siRNA)  %>%
  summarise(meanCt = mean(Ct))

# Join the dataset with the housekeeping gene information from the same treatment group
# Drop Ct and housekeeping meanCt columns
FibrosisArraydCt <- left_join(FibrosisArrayCt, SummaryHouseKeeping, by = c("Time", "n", "Agonist", "siRNA")) %>%
  mutate(dCt = Ct - meanCt) %>%
  select("Gene", "Time", "Agonist", "siRNA", "n", "dCt")

# Find average dCt for the control condition
SummaryControl <- FibrosisArraydCt %>%
  filter(siRNA == "siCtrl",
         Agonist == "DMEM") %>%
  select("Gene", "Time", "Agonist", "n", "dCt")

# Join the control dCt value with the rest of the dataset
FibrosisArrayddCt <- left_join(FibrosisArraydCt, SummaryControl, by = c("Gene", "Time", "n")) %>%
  rename(Agonist = Agonist.x) %>%
  # Calculate the fold change using the the 2^(-ΔΔCt)
  mutate(ddCt = dCt.x - dCt.y, 
         FC = 2^-ddCt,
         ddCt = -ddCt) %>%
  select("Gene", "Time", "Agonist", "siRNA", "n", "ddCt", "FC") %>%
  group_by(Gene, Time, Agonist, siRNA) %>%
  # Find mean ddCt and fold change (FC) for each gene across the biological replicates within each treatment
  summarise(meanddCt = mean(ddCt),
            meanFC = mean(FC))

# Creaete boxplot of mean ddCt for all genes within a treatment
ggplot(FibrosisArrayddCt, aes(x = Agonist, y = meanddCt, color = siRNA)) +
  geom_boxplot() +
  facet_wrap(~Time) +
  labs(y = "log2(Fold Change) \n (normalized to control)") +
  scale_color_manual(values=c("#000000", "#4387b5", "#bf3c3c")) +
  theme_pubclean() + 
  theme(legend.position="bottom",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line.y.left = element_line()) +
  geom_hline(yintercept = 0, colour = "black", size = 0.5, linetype = "dashed")

# Perform select t-test with Bonferonni correction for every gene with the dCt values
# Modified t-test with pooled variance and select comparisons
t_test_select <- function(y, test_list) {
  x <- y$dCt
  g <- y$Treatment
  g <- factor(g)
  xbar <- tapply(x, g, mean, na.rm = TRUE)
  s <- tapply(x, g, sd, na.rm = TRUE)
  n <- tapply(!is.na(x), g, sum)
  degf <- n - 1
  total.degf <- sum(degf)
  pooled.sd <- sqrt(sum(s^2 * degf)/total.degf)
  p.values <- list()
  test_list <- test_list
  if (is.list(test_list)) {
  for(var in test_list) {
    dif <- xbar[var[1]] - xbar[var[2]] 
    se.dif <- pooled.sd * sqrt(1/n[var[1]] + 1/n[var[2]])
    t.val <- dif/se.dif
    p.values[[length(p.values) + 1]] <- 2 * pt(-abs(t.val), total.degf)
  }
  } else {
    print('Test list is not properly formatted')
  }
  names(p.values) <- test_list
  p.values <- enframe(p.values) 
  p.values <- unnest(p.values, value)
  adj <- p.adjust(p.values$value, method = "bonferroni")
  cbind(p.values, adj)
}

# Split dataset by a list with a dataset for each gene
FibrosisArraydCt_test <- FibrosisArraydCt %>%
  # Join the aspects of the treatments into one column
  unite(Treatment, Agonist, siRNA, sep = "/") %>%
  split(., f = .$Gene) %>%
  # Perform a t-test between treatment groups of choosing
  # provide a vector comprised of the two names of the two treaetments
  map_dfr(t_test_select, 
          .id = 'names',
          test_list = list(c("treatment1", "treatment2"),
                           c("treatment2", "treatment3")))




