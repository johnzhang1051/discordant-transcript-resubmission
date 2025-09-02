any_two<- read.csv("Data/any_two.csv")
unique_promoter <- read.csv("Data/unqiue_promoter.csv")
final <- merge(any_two, unique_promoter, by ="transcript_id")

any_one <- read.csv("Data/transcript_unique_with_any_peak.csv")
final_any_one <- final <- merge(any_one, unique_promoter, by ="transcript_id")