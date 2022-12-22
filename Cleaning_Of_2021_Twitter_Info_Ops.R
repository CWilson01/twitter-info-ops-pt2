######Package load########

# Loads packages using the groundhog package. Groundhog enables reproducible     
# analysis by recording the date the packages were used. It then downloads those 
# exact package version upon later analysis, ensuring scripts run as intended.

library("groundhog")
pkgs <- c("tidyverse", "skimr")
groundhog.library(pkgs, "2022-07-11")

# Disabling scientific notation to improve readability of certain variables.
# options(scipen=0, digits=7) # To return to default setting
options(scipen = 999)

# Setting seed for results replication
set.seed(12345)

######Data Load##########

# Load the three datasets of removed tweets acquired from Twitter. Feb 2021 removal
# of Russian IRA accounts, Feb 2021 removal of Iranian accounts, and Dec 2021
# removal of PRC Xinjiang-focused accounts.

russia_2021 <- read_csv("Russia 2021/Feb 2021/hashed_2020_12_IRA_202012_IRA_202012_tweets_csv_hashed.csv")
iran_2021<- read_csv("Iran 2021/hashed_2020_12_iran_202012_iran_202012_tweets_csv_hashed.csv")
china_2021 <- read_csv("China 2021/Peoples_republic_of_China_Xinjiang_DEC2021_tweets/CNHU_0621_tweets_csv_hashed_2021.csv")

## Examine raw data
skim_without_charts(russia_2021)
skim_without_charts(iran_2021)
skim_without_charts(china_2021)

######Data Cleaning######

# Transform tweet and account identifiers from numeric to character.
# Intended to ensure that these are used for identification purposes rather than 
# descriptive statistics.
russia_2021 <- russia_2021 %>%
  mutate_at(
    c('tweetid',
    'in_reply_to_userid',
    'in_reply_to_tweetid',
    'quoted_tweet_tweetid',
    'retweet_userid',
    'retweet_tweetid',
    'user_mentions'),
    as.character)

iran_2021 <- iran_2021 %>%
  mutate_at(
    c('tweetid',
      'in_reply_to_userid',
      'in_reply_to_tweetid',
      'quoted_tweet_tweetid',
      'retweet_userid',
      'retweet_tweetid',
      'user_mentions'),
    as.character)

china_2021 <- china_2021 %>%
  mutate_at(
    c('tweetid',
      'in_reply_to_userid',
      'in_reply_to_tweetid',
      'quoted_tweet_tweetid',
      'retweet_userid',
      'retweet_tweetid',
      'user_mentions'),
    as.character)

# Removing case sensitivity from hashtags variable to ensure consistent counting
# of hashtag usage
russia_2021$hashtags <- tolower(russia_2021$hashtags)
iran_2021$hashtags <- tolower(iran_2021$hashtags)
china_2021$hashtags <- tolower(china_2021$hashtags)

# Removing brackets and single quotes surrounding strings in two columns
russia_2021$user_mentions <- gsub("\\[|\\]", "", russia_2021$user_mentions)
russia_2021$user_mentions <- gsub("'", "", russia_2021$user_mentions)
russia_2021$hashtags <- gsub("\\[|\\]", "", russia_2021$hashtags)
russia_2021$hashtags <- gsub("'", "", russia_2021$hashtags)

iran_2021$user_mentions <- gsub("\\[|\\]", "", iran_2021$user_mentions)
iran_2021$user_mentions <- gsub("'", "", iran_2021$user_mentions)
iran_2021$hashtags <- gsub("\\[|\\]", "", iran_2021$hashtags)
iran_2021$hashtags <- gsub("'", "", iran_2021$hashtags)

china_2021$user_mentions <- gsub("\\[|\\]", "", china_2021$user_mentions)
china_2021$user_mentions <- gsub("'", "", china_2021$user_mentions)
china_2021$hashtags <- gsub("\\[|\\]", "", china_2021$hashtags)
china_2021$hashtags <- gsub("'", "", china_2021$hashtags)

# Adding NAs to missing observations in the above columns
russia_2021$user_mentions <- na_if(russia_2021$user_mentions, "")
russia_2021$hashtags <- na_if(russia_2021$hashtags, "")

iran_2021$user_mentions <- na_if(iran_2021$user_mentions, "")
iran_2021$hashtags <- na_if(iran_2021$hashtags, "")

china_2021$user_mentions <- na_if(china_2021$user_mentions, "")
china_2021$hashtags <- na_if(china_2021$hashtags, "")

## Examine cleaned data
skim_without_charts(russia_2021)
skim_without_charts(iran_2021)
skim_without_charts(china_2021)
# no missing values were noted in any variables of interest, cleaning is completed

######Writing Data########
# Write cleaned datasets to new files to analysis
write_csv(russia_2021, "russia_2021_cleaned.csv")

write_csv(iran_2021, "iran_2021_cleaned.csv")

write_csv(china_2021, "china_2021_cleaned.csv")
