source("https://raw.githubusercontent.com/walmes/wzRfun/master/R/pairwise.R") # Cria os contrastes 2 x2. Uso a função acp()
## --------------------------------------------------------------------------
# Interação entre 2 categóricas
## --------------------------------------------------------------------------

# fixed <- "estac * season"
# effect <- c("estac", "season")
# Matriz para interação 2x2 - Categóricas
mat <- function(fixed, effect, data){
  data$y <- runif(nrow(data))
  form <- as.formula(paste("y", " ~ ", fixed))
  m0 <- lm(form, data = data)
  K <- doBy::LE_matrix(m0, effect = effect)
  grid <- attr(K, "grid")
  grid <- grid[, 1:2]
  K2 <- by(K, grid[, 1], as.matrix)                       # O 1º efeito é fixado; e o 2º efeito é comparado dentro das linhas
  K2 <- lapply(K2, FUN = "rownames<-", unique(grid[, 2])) # Somente coloca os nomes nas linhas
  return(K2)
}
# k1 <- mat(fixed, effect, m1$data)
# Interação 2x2
glht <- function(mat, model, resp = 1, alpha = 0.05, transf = F){
  coef.mat <- coef(model, type = "beta")      # Obtem TODOS os coeficientes fixos do modelo
  coefs <- coef.mat[coef.mat$Response==resp,] # Obtem da variável resposta em questão
  positions <- as.numeric(row.names(coef.mat[coef.mat$Response==resp, ])) #Posições das variáveis
  vcovs <- as.matrix(vcov(model))[positions , positions] # Matriz de variância e covariância para a variável resposta
  
  # Medidas Pontuais para cada parâmetros
  mu <- mat%*%coefs$Estimates      # Média na Escala do Preditor Linear
  var <- mat%*%vcovs%*%t(mat)      # Variância na Escala do Preditor Linear
  e.p <- sqrt(diag(var))           # Erro Padrão dos coeficientes na escala do Preditor Linear
  quantil <- 1 - (alpha/2)
  # Calcula as médias ajustadas de cada parâmetro
  mu.int <- data.frame(Media = mu,
                       ErroPadrao   = e.p,
                       ic.inf = mu-qnorm(quantil)*e.p,  # IC via distribuição Z
                       ic.sup = mu+qnorm(quantil)*e.p  # IC via distribuição Z
  )
  
  # Calcula os Contrastes
  ctr <- apc(mat) # Essa função do prof. Walmes calcula todas as comparacoes 
  d <- ctr%*%coefs$Estimates   # Calcula a diferença em média na escala do preditor
  var <- ctr%*%vcovs%*%t(ctr)  # Calcula a variância da diferença em média na escala do preditor
  sd <- sqrt(diag(var))        # Calcula o erro padrão da diferença em média na escala do preditor
  
  # Formata, calcular medidas adicionais
  z.value <- abs(d/sd)
  p.values <- (2*pnorm(z.value, lower.tail = F)) # Bilateral
  p.values <- ifelse(p.values>=1,1,p.values)
  p.values.adj <- p.adjust(p.values, method = "fdr")
  p.values.adj <- ifelse(p.values.adj>=1,1,p.values.adj)
  # Intervalos de confiança
  IC.Inf <- d - qnorm(quantil)*sd # Intervalo de confiança sem correção
  IC.Sup <- d + qnorm(quantil)*sd
  resumo <- data.frame(Estimativa = d,
                       `Erro Padrão` = sd,
                       Valor.Z = z.value,
                       Valorp = p.values,
                       Valorp.adj = round(p.values.adj,6),
                       IC.Inf = IC.Inf,
                       IC.Sup = IC.Sup)
  
  if (transf){
    mu.int <- within(mu.int, {
      Media <- exp(Media)
      ic.inf <- exp(ic.inf)
      ic.sup <- exp(ic.sup)
    })
    mu.int <- mu.int[, -2] # Retira o Erro padrão pois não está na escala do log
    resumo <- within(resumo, {
      Estimativa <- exp(Estimativa)
      IC.Inf <- exp(IC.Inf)
      IC.Sup <- exp(IC.Sup)
    })
    resumo <- resumo[, -2]
  }
  print("A função realiza as comparações múltiplas da 2a variável fixando um nível da primeira variável. O nível fixado da primeira variável corresponde ao índice que está na matriz K")
  print("Correção de fdr para o valor-p; Intervalo de Confiança nominal de 95%")
  list(glht = resumo,
       means = mu.int)
}
# glht(k1[[5]], m1, resp = 3, transf = F)

## --------------------------------------------------------------------------
# Interação entre categórica e numérica
## --------------------------------------------------------------------------
# model <- m1
# resp <- 6
# data <- dat
# fixed <- "estac*ano"
intcf <- function(model, resp, data, fixed, alpha = 0.05){
  coef.mat <- coef(model, type = "beta") #Obtem TODOS os coeficientes do modelo
  positions <- as.numeric(row.names(coef.mat[coef.mat$Response == resp,])) #Posições das variáveis
  coefs <- coef.mat[positions, ]
  vcovs <- as.matrix(vcov(model))[positions , positions] # Matriz de variância e covariância para a variável resposta
  # Versão 1 - Modelo
  form <- as.formula(paste("~ ", fixed))
  Xp <- model.matrix(form, data)
  se2 <- unname(rowSums((Xp %*% vcovs) * Xp))
  quantil <- 1 - alpha/2
  qt <- qnorm(quantil)
  
  data$pred <- as.numeric(Xp%*%coefs$Estimates)
  data$CI.inf <- (data$pred + -qt*sqrt(se2))
  data$CI.sup <- (data$pred + qt*sqrt(se2))
  
  data$pred.exp <- exp(data$pred)
  data$CI.inf.exp <- exp(data$pred + -qt*sqrt(se2))
  data$CI.sup.exp <- exp(data$pred + qt*sqrt(se2))
  
  return(data)
}
# dat <- intcf(m1, 6, dat, "estac * ano")


## --------------------------------------------------------------------------
# 1 covariável categórica
## --------------------------------------------------------------------------
mat2 <- function(fixed, data){
  data$y <- runif(nrow(data))
  form <- as.formula(paste("y", " ~ ", fixed))
  m0 <- lm(form, data = data)
  K <- model.matrix(m0)
  K <- unique(K)
  rownames(K) <- levels(data[, fixed]) # Somente coloca os nomes nas linhas
  return(K)
}
# Inverse logit function for binomial family
inv.logit <- function(x){
  1/(1+exp(-x))
}

# Formata o resultado do glht.s
formatar <- function(a, resp){
  print(kable(a[[1]], 
        caption = paste0("Proporção da proteína ", resp, " ajustada pelo modelo para cada tipo de corte"), 
        align = "c"))

  print(kable(a[[2]], 
        caption = paste0("Razão de Chances entre os tipos de corte para a proteína ", resp),
        align = "c"))
}
# fixed <- "estac" # Variável que você quer comparar
# mat <- mat2(fixed, dat) # Matriz

# resp <- 1        # Ordem da resposta do mcglm
# order <- 2:5     # Ordem dos coeficientes "fixed"

# Comparações para 1 único fator
glht.s <- function(model, resp, covariate, data, order, transf = F, alpha = 0.05){
  mat <- mat2(covariate, data)
  coef.mat <- coef(model, type = "beta")      # Obtem TODOS os coeficientes fixos do m1o
  coefs <- coef.mat[coef.mat$Response==resp,] # Obtem da variável resposta em questão
  positions <- as.numeric(row.names(coef.mat[coef.mat$Response==resp,])) #Posições das variáveis
  vcovs <- as.matrix(vcov(model))[positions , positions] # Matriz de variância e covariância para a variável resposta
  coefs <- coefs[c(1, order), ]
  vcovs <- vcovs[c(1, order), c(1, order)]
  # Medidas Pontuais para cada parâmetros
  mu <- mat%*%coefs$Estimates      # Média na Escala do Preditor Linear
  var <- mat%*%vcovs%*%t(mat)      # Variância na Escala do Preditor Linear
  e.p <- sqrt(diag(var))           # Erro Padrão dos coeficientes na escala do Preditor Linear
  quantil <- 1 - (alpha/2)
  # Calcula as médias ajustadas de cada parâmetro
  mu.int <- data.frame(Media = mu,
                       ErroPadrao   = e.p,
                       ic.inf = mu-qnorm(quantil)*e.p,  # IC via distribuição Z
                       ic.sup = mu+qnorm(quantil)*e.p)  # IC via distribuição Z
  ctr <- apc(mat) # Essa função do prof. Walmes calcula todas as comparacoes 
  d <- ctr%*%coefs$Estimates   # Calcula a diferença em média na escala do preditor
  var <- ctr%*%vcovs%*%t(ctr)  # Calcula a variância da diferença em média na escala do preditor
  sd <- sqrt(diag(var))        # Calcula o erro padrão da diferença em média na escala do preditor
  
  # Formata, calcular medidas adicionais
  z.value <- abs(d/sd)
  p.values <- (2*pnorm(z.value, lower.tail = F)) # Bilateral
  p.values <- ifelse(p.values>=1,1,p.values)
  p.values.adj <- p.adjust(p.values, method = "fdr")
  p.values.adj <- ifelse(p.values.adj>=1,1,p.values.adj)
  # Intervalos de confiança
  IC.Inf <- d - qnorm(quantil)*sd # Intervalo de confiança sem correção
  IC.Sup <- d + qnorm(quantil)*sd
  resumo <- data.frame(Estimativa = d,
                       `Erro Padrão` = sd,
                       Valor.Z = z.value,
                       Valorp = p.values,
                       Valorp.adj = round(p.values.adj,6),
                       IC.Inf = IC.Inf,
                       IC.Sup = IC.Sup)
  
  if (transf == "exp"){
    # Transformation of the estimates
    mu.int <- within(mu.int, {
      Media <- exp(Media)
      ic.inf <- exp(ic.inf)
      ic.sup <- exp(ic.sup)
    })
    mu.int <- mu.int[, -2] # Retira o Erro padrão pois não está na escala do log
    # Transformation of the contrasts
    resumo <- within(resumo, {
      Estimativa <- exp(Estimativa)
      IC.Inf <- exp(IC.Inf)
      IC.Sup <- exp(IC.Sup)
    })
    resumo <- resumo[, -2]
  } else if (transf == "inv.logit"){ # For proportion; 0/1 data
    # Transformation of the estimates
    mu.int <- within(mu.int, {
      Media <- inv.logit(Media)*100
      ic.inf <- inv.logit(ic.inf)*100
      ic.sup <- inv.logit(ic.sup)*100
    })
    mu.int <- mu.int[, -2] # Retira o Erro padrão pois não está na escala da prob
    names(mu.int) <- c("Média (%)", "IC.Inf 95% (%)", "IC.Sup 95% (%)")
    # Transformation of the contrasts
    resumo <- within(resumo, {
      Estimativa <- inv.logit(Estimativa)
      IC.Inf <- inv.logit(IC.Inf)
      IC.Sup <- inv.logit(IC.Sup)
    })
    resumo <- resumo[, -2] # Retira o Erro padrão pois não está na escala de prob
    names(resumo) <- c("Razão de Chances", "Z", "Valor-p", "Valor-p-FDR", "IC.Inf 95%", "IC.Sup 95%")
  }
  list(mu.int,
       resumo)
}
# glht.s(m1, resp, mat, 2:5)
