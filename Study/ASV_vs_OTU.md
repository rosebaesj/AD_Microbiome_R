# ASV vs OTU

by Sunjeong Bae \
date 2022 May 13 \
2022 May 12 일 회의 후 ASV 와 OTU에 대한 충분한 이해가 필요하다 생각되어 정리하는 자료 \
참고: https://www.zymoresearch.com/blogs/blog/microbiome-informatics-otu-vs-asv \

## 과거부터 사용되어오는 방법인 OTU

OTU도 두가지 방법이 있다.

### 1. Reference-free OTU Clustering

de novo. \
즉 machine learning 방법 등을 활용하여 알아서 끼리끼리 묶는것이라 보면 될듯하다. \
우리가 사용할 일은 없는 듯 하다. \

### 2. Reference-based OTU Clustering

이게 이제 보통 사용하는 OTU Clustering이다 \
database가 있으면 거기에서 비슷한것끼리 묶어 넣는다. \

## Amplicon Sequence Variant (ASV)

진짜 하나의 read가 하나의 feature\
diversity를 더 정확히 분석할 수 있다.\
ASV 분석 중 가장 sensitive 한것이 DADA2\

### Chimera 문제
Because an ASV is an exact sequence, a chimeric ASV can be expected to be the exact or near-exact child of two more prevalent exact parent sequences in the same sample, with one parent contributing the left side and one parent contributing the right side of the chimera. These properties can enable the identification of chimeras by using the alignment of less prevalent ASVs to more prevalent ASVs from the same sample8.

### 아니 그래서 ASV 에 taxonomy 정보 입힐때는 걍 되는겨? 뭔겨.....
