data_tofix <- read.csv("Intern_intermed.csv", header = T, stringsAsFactors = F) %>%
  mutate(other = as.numeric(other), sculp = as.numeric(sculp), gun = as.numeric(gun)) %>%
  filter(name != "") 
data_tofix[is.na(data_tofix)] <- 0

new_data <- data_tofix %>%
  group_by(date, name) %>%
  summarize(gun = sum(gun), sculp = sum(sculp), other = sum(other))

write.csv(new_data, "new_data.csv", row.names = F)

