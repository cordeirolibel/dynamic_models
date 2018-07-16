# III Escola de Matemática Aplicada
# MODELOS DINÂMICOS BAYESIANOS

Projeto desenvolvido na [III Escola de Matemática Aplicada](http://www.cemeai.icmc.usp.br/component/jem/event/65-3-escola-de-matematica-aplicada) do curso [Modelos Dinâmicos e Aplicaçes](http://www.cemeai.icmc.usp.br/component/k2/item/681-modelos-dinamicos-e-aplicacoes) ministrado pelo professor [Helio S. Migon](http://www.dme.ufrj.br/migon/) (UFRJ). <br>
__Resumo:__ A Modelagem dinâmica e previsão  Bayesiana  de séries temporais é uma das mais importantes áreas surgidas no final do século passado.
Tem como ponto de partida os modelos de regressão. A inclusão de uma evolução temporal dos parâmetros da regressão dão origem a uma ampla classe de modelos incluindo vários temas da estatística moderna. (ver Migon, H.S et al,Handbook ofStatistics Vol 25, Cap. 19, 553-588, 2006. e também Pole West e Harrison, 1994).

__Ementa:__
1. Introdução aos modelos dinâmicos.<br>
2. Modelos dinâmicos lineares.<br>
3. Aspectos teóricos de modelos dinâmicos de séries temporais.<br>
4. Modelos polinomiais, sazonais e  de regressão.<br>
5. Monitoramento e intervenção.<br>
6. Tópicos especiais: MLG dinâmicos,  Modelos Hierárquicos dinâmicos.<br>
7. Exemplos: Chuva-vazão, propaganda em TV, tábuas de mortalidade, modelos de credibilidade dinâmico etc.

-----------------------------
# Projeto 
&emsp;O objetivo do projeto consiste em consolidar os conhecimentos do curso em uma aplicação. 
  Cada grupo de alunos escolhia uma base dados para aplicar os modelos.

__Grupo:__
Emanuelle A. Paixão (UFLA) - 
Gustavo C. Libel (UTFPR) - 
Junior R. Ribeiro (ICMC)

__Descrição dos dados utilizados:__
Os dados utilizados para análise correspondem a uma série mensal de Janeiro/2004 a Julho/2018, descrevendo a porcentagem de pesquisa no Google do termo “USP” acessados no Brasil, com dados extraídos do Google Trends (Figura 1 - (a) e (b)). Os dados a serem analisados serão os constantes na Figura 1 (a), sendo a Figura 1 (b) apenas ilustrativa.<br>
&emsp; OBSERVAÇÃO: Consideramos que os dados apresentam uma distribuição normal.

|<img src="imgs/usp.PNG">|
|:--:| 
|*(a) Percentual de acesso do termo “USP” de janeiro de 2004 a julho de 2018.*|
|<img src="imgs/5anos.PNG">|
|*(b) Percentual de acesso do termo “USP” de julho de 2013 a julho de 2018.*
|Figura 1: Percentual de pesquisas do termo “USP” no Google no Brasil.|

&emsp; Na Figura 1 constam dados extraídos da plataforma do Google Trends, contabilizando a quantidade de pesquisas do termo “USP” no Google Web ao longo dos anos desde janeiro de 2004 até julho de 2018. Os dados capturados estão normalizados pelo maior valor (frequência relativa).<br>
&emsp; Buscando uma maior riqueza de detalhes, colocamos na Figura 1 (b) um zoom in sobre o gráfico, ou seja, capturando os dados a partir de julho de 2013 até julho de 2018. Nesse gráfico, pode-se perceber que no mês de dezembro há uma baixa frequência de pesquisas do termo (vales acentuados), possivelmente causada por conta do período de recesso escolar. Em contrapartida, nos meses seguintes, janeiro e fevereiro, nota-se um aumento expressivo no interesse do termo “USP”, causado provavelmente pela publicação dos resultados dos processos seletivos e vestibulares e em seguida a matrícula de ingressantes na universidade.<br>
&emsp; De uma forma geral, pode-se perceber, um decaimento acentuado no número de pesquisas do termo “USP” ao longo dos anos, com alguns picos e vales, o que pode caracterizar a existência de uma possível sazonalidade. O decaimento talvez sugere que, quando a internet começou a se popularizar, as pessoas tinham grande curiosidade em conhecer a universidade, curiosidade essa que foi diminuindo conforme a divulgação e difusão de oportunidades pelo meio eletrônico, por ser uma realidade mais rotineira. Percebe-se que os maiores valores de pico ocorrem geralmente nos meses de fevereiro e agosto, o que pode ser explicado pelas matrículas da universidade, que são semestrais, e então reflete em um maior número de acessos. Os menores valores (vales) são observados nos finais de semestre, que são em julho e dezembro, período de recesso escolar. Nos gráficos mais adiante, essa sazonalidade será mostrada com mais detalhes.<br>
&emsp; Na Figura 2 apresentamos o ajuste do Modelo Linear de Primeira Ordem e na Figura 3, o Modelo Linear de Segunda Ordem. Podemos observar a existência de diversos pontos fora do intervalo de confiança e também que a curva de predição não descreve totalmente a variação dos dados. Em vermelho consta a previsão em 1 passo, em azul as faixas de confiança e os pontos são os dados.<br>

|<img src="imgs/primeira_ordem.jpeg">|
|:--:| 
|*Figura 2: Ajuste do Modelo Dinâmico de Primeira Ordem aos dados analisados, considerando m0 = 0, C0 = 1000, V = 400 e W = 4.*|

|<img src="imgs/segunda_ordem.jpeg">|
|:--:| 
|*Figura 3: Ajuste do Modelo Dinâmico de Segunda Ordem aos dados analisados, considerando m0 = c(0,0), C0 = diag(1000,n,n) e W = diag(c(0.02, 0.01)) (a) com o respectivo fator de crescimento (b).*|

&emsp; Como pode-se perceber ao comparar as Figuras 2 e 3, o ajuste de um Modelo Linear de Segunda Ordem pode descrever melhor as variações dos dados, pois menos pontos ficaram fora da faixa de confiança.
&emsp; A sazonalidade dos dados pode ser modelada usando séries harmônicas truncadas. Na Figura 4, é apresentada uma modelagem considerando uma série harmônica truncada em apenas 1 harmônico.

|<img src="imgs/1_harm.png">|
|:--:| 
|*Figura 4: Ajuste do Modelo Dinâmico com 1 Harmônico aos dados analisados, considerando w=2π/12, D1=matrix(1/sqrt(.8),2), D2=matrix(1/sqrt(.95),2), D=bdiag(D1,D3), G1=matrix(c(1,0,1,1),2), G2=matrix(c(cos(w),-sin(w),sin(w),cos(w)),2), G=bdiag(G1, G2),m0=rep(0,n), C0=diag(100,n), n0=10, d0=1 e F=cbind(matrix(c(rep(1,T),rep(0,T),rep(1,T),rep(0,T)),nrow=n,byrow=TRUE),c(1,0,1,0)).*|

&emsp; O efeito do harmônico sobre os dados é mostrado na Figura 5. Nota-se uma suavização ao longo do tempo por conta de que a variância dos dados foi diminuindo.

|<img src="imgs/fat_harm.png">|
|:--:| 
|*Figura 5: Efeito com 1 harmônico sobre os dados.*|

__REFERÊNCIAS:__

____ Plataforma Google Trends. https://trends.google.com.br/trends. Acesso em 14 de julho de 2018. <br>
[MAROTTA, Raíra.](https://github.com/rairamarotta) Algoritmo usado durante o curso.
