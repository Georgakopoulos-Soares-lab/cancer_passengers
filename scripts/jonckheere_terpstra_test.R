library(clinfun)

inp_file = "results/median_survival_time_pancancer_groups.csv"
inp_data = read.csv(inp_file, header = TRUE)
head(inp_data)

group = factor(inp_data$group, levels = c("Has no mutations", "Low", "Medium", "High"), ordered = TRUE)
values = inp_data$median_survival
print(head(group))

jonckheere.test(values, group)