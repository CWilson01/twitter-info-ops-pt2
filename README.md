# Under the Radar: Analyzing Recent Twitter Information Operations to Improve Detection and Removal of Malicious Actors, Part 2
R code, markdown document, and Gephi files used to create a social network analysis project examining three information operations removed from Twitter in 2021.

The full written report of my findings can be found on my website, https://wonksecurity.com. Here is the direct link to the pdf of the report: https://wonksecurity.com/wp-content/uploads/2022/12/Twitter_info_op_report_v2.pdf. Additionally, you can also find part 1 of this project [here](https://github.com/CWilson01/twitter-info-ops-pt1).

## Abstract
This report builds upon the work done in part one of this series by examining the network structure of three information operations (IOs) that were removed from Twitter in 2021.  The analysis that follows uses social network analysis (SNA) to explore the structure, key network statistics, and measures of centrality for network graphs created from Twitter mentions. Five data sets feature in this analysis, three IO networks and two control networks. Data for the three IO networks came from Twitter’s Transparency Center and contained tweets from a Russian, Chinese, and Iranian IO, respectively. In addition, two COVID-19 tweet data sets from Kaggle served as the controls. This project seeks to determine if it is possible to make cross-network comparisons that could enhance the early detection of IOs on social media platforms like Twitter. The analysis found that while each network was structurally unique, the key SNA statistics failed statistical significance testing when checking for differences between the IO group and control group. This may be the result of a small network sample size (n=5). However, this study also found that measures of centrality had statistically significant differences between the IO group and the control group. This suggests that measures of centrality, particularly eigenvector centrality and Pagerank, could be useful metrics for differentiating IOs from legitimate Twitter conversations.

## How to use this project:

- Go to [Twitter's Transparency Center](https://transparency.twitter.com/en/reports/moderation-research.html) and acquire the original datasets used in this project. Given the rapidly evolving policies and changes at Twitter currently, I will not be making the raw data available while it remains easily accessible via Twitter's own service.

- Acquire the "People’s Republic of China - Xinjiang (December 2021) - 2048 Accounts" dataset released in December 2021, particularly the Tweet Information file.

- Acquire the Tweet Information file from the "Iran (February 2021) - 238 Accounts" dataset released in February 2021.

- Acquire the Tweet Information file from the "Russia IRA (February 2021) - 31 Accounts" dataset also releaed in February 2021.

- Then, go to Kaggle user [Arunava Kumar Chakraborty's COVID-19 dataset](https://www.kaggle.com/datasets/arunavakrchakraborty/covid19-twitter-dataset) page and download the 61MB file that contains both COVID-19 datasets that were used as controls.

- Download the R scripts and Quarto markdown file, and Gephi files provided here to replicate, build upon, or fork the cleaning and analysis process that I used.

These datasets provided a massive wealth of possible information, and it is my goal to conduct additional analysis in future projects. I list out several potential ideas of how the data could be further analyzed in part one of this project, located [here](https://github.com/CWilson01/twitter-info-ops-pt1). If you try one of these or have additional ideas, feel free to drop me a line at cody [@] wonksecurity [.] com. I'd love to hear what you think.

## License
The code, markdown, and Gephi files in this project are released under a GPL-3.0 license. The data from Twitter's Transparency Center is bound by its terms of use, found [here](https://developer.twitter.com/en/developer-terms). The COVID-19 datasets made available by A. K. Chakraborty are available on Kaggle under a [CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0/) license. Full credit for this dataset goes to: Chakraborty A.K., Das S., Kolya A.K. (2021) Sentiment Analysis of Covid-19 Tweets Using Evolutionary Classification-Based LSTM Model. In: Pan I., Mukherjee A., Piuri V. (eds) Proceedings of Research and Applications in Artificial Intelligence. Advances in Intelligent Systems and Computing, vol 1355. Springer, Singapore. https://doi.org/10.1007/978-981-16-1543-6_7.
