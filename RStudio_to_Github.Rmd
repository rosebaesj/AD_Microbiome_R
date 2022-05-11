---
title: "RStudio_to_Github"
output: html_document
date: '2022-05-11'
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## R Markdown

This is an R Markdown document. Markdown is a simple formatting syntax for authoring HTML, PDF, and MS Word documents. For more details on using R Markdown see <http://rmarkdown.rstudio.com>.

When you click the **Knit** button a document will be generated that includes both content as well as the output of any embedded R code chunks within the document. You can embed an R code chunk like this:

```{r cars}
summary(cars)
```

## Including Plots

You can also embed plots, for example:

```{r pressure, echo=FALSE}
plot(pressure)
```

Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.




이후 github token을 아래 콘솔에서 입력해주면 됨.

편집이 끝난 것들은 commit을 통해 저장을 해야함
이후 실제 원격 저장소에 넣기 위해서는 초록색 화살표로 push!
오른쪽 위 Git에서 push화살표를 눌러주면 github로 동기화 완료!!

ㅊCommit과 push의 차이점
https://rateye.tistory.com/1916
그리고나서 오른쪽 Git에서 push화살표를 눌러주면 github로 동기화 완료!!

Git  ! [rejected]        HEAD -> main (non-fast-forward) 의 오류가 나는 ㅕㄱㅇ우
Terminal에서 
$ git push origin +main
진행하면 main과 안맞는 부분이 있더라도 강제 덮어쓰기....
문제가 될 것 같기는 함.... 해결 방법 모르겠음
그때그때 pull 잘하고 push 하자.

Pull에서 오류가 있는 것 같다.
hint:   git config pull.rebase false  # merge (the default strategy)
hint:   git config pull.rebase true   # rebase
hint:   git config pull.ff only       # fast-forward only




