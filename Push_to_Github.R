#################################################
#*** Push to Github ***#
#################################################
#* 이 R code는 무시하고 RStudio_to_Github.Rmd를 확인하기 바람

install.packages("gitcreds")
library(gitcreds)
gitcreds_set()
#*이후 github token을 아래 콘솔에서 입력해주면 됨.
#*
#*편집이 끝난 것들은 commit을 통해 저장을 해야함
#*이후 실제 원격 저장소에 넣기 위해서는 초록색 화살표로 push!
#*오른쪽 위 Git에서 push화살표를 눌러주면 github로 동기화 완료!!
#*
#*ㅊCommit과 push의 차이점
#*https://rateye.tistory.com/1916
#*그리고나서 오른쪽 Git에서 push화살표를 눌러주면 github로 동기화 완료!!
#*
#*Git  ! [rejected]        HEAD -> main (non-fast-forward) 의 오류가 나는 ㅕㄱㅇ우
#*Terminal에서 
#*$ git push origin +main
#*진행하면 main과 안맞는 부분이 있더라도 강제 덮어쓰기....
#*문제가 될 것 같기는 함.... 해결 방법 모르겠음
#*그때그때 pull 잘하고 push 하자.
#*
#*Pull에서 오류가 있는 것 같다.
#*hint:   git config pull.rebase false  # merge (the default strategy)
#*hint:   git config pull.rebase true   # rebase
#*hint:   git config pull.ff only       # fast-forward only
#*
#*
#*
#*
