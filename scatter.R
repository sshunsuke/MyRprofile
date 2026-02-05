# scatter1
# simple scatter diagram

x <- c(1, 2, 3, 4, 5)
y <- c(5, 7, 3, 8, 2)

#pdf(file="./scatter1_(simple).pdf")
png(file="./scatter1_(simple).png")
plot(x, y, main="Scatter 1 (Simple)", xlab="X Label", ylab="Y Label", pch=15)
dev.off()


# - - - - - - - - - - - - - - - - - - - -
# scatter2
# multi-data

scatter2 <- function() {
  x <- c(0, 6, 9, 10, 20)
  y1 <- c(5, 7, 3, 8, 2)
  y2 <- c(15, 11, 13, 19, 7)
  y3 <- c(24, 21, 27, 1, 15)

  xRange <- c(0, 20)
  yRange <- c(0, 30)

  # 空のグラフ
  plot(xRange, yRange, type = 'n', main="Scatter 2 (multi-data)", xlab="x", ylab="y")

  # 格子の線を引く場合  (大雑把で良ければ、grid(5,5) みたいな感じでよい)
  abline(v = 0:4*5, h = 0:6*5, col = "gray", lty = 2, lwd = 0.7)

  # 実データをプロット
  lines(x, y1, col="black", lty=1, lwd=4)
  lines(x, y2, col="red",   lty=2, lwd=3)
  lines(x, y3, col="blue",  lty=3, lwd=2)

  # 凡例を入れる
  legend(15, 25, c("A", "B", "C"), lty=c(1, 2, 3), lwd=c(4, 3, 2),
         col=c("black", "red", "blue"), bg="white")

  # 説明テキスト
  text(15, 26, "Text", adj=0, font=2, cex=0.9, col="#0080FF" )
}

#pdf(file="./scatter2_(multi-data).pdf")
png(file="./scatter2_(multi-data).png")
scatter2()
dev.off()


# - - - - - - - - - - - - - - - - - - - -

# scatter3
# time-line

scatter3 <- function() {
  x <- as.POSIXct( c("2012-02-12 10:0:0", "2012-02-12 12:30:0", "2012-02-12 15:0:0", "2012-02-12 20:0:0") )
  y <- c(2, 5, 7, 3)
  xRange <- c( x[1], tail(x, 1) )
  yRange <- c( 0, 10 )

  plot(xRange, yRange, type='n', xaxt="n", main="Scatter 3 (axis-time)", xlab="Time", ylab="Data")
  axis.POSIXct(1, at=seq(xRange[1], xRange[2], by="2 hour"), format="%H:%M")
  abline(v = seq(x[1], tail(x, 1), by="2 hour"), col = "#dddddd", lty = 1, lwd = 0.7)
  lines(x, y)
  points(x, y, pch=19)
}

#pdf(file="./scatter3_(axis-time).pdf")
png(file="./scatter3_(axis-time)).png")
scatter3()
dev.off()


# - - - - - - - - - - - - - - - - - - - -
# scatter4
# 2軸 (Y)

scatter4 <- function() {
  x <- c(0, 10, 20, 30, 50)
  y1 <- c(0, 1, 2, 4, 8)
  y2 <- c(900, 450, 225, 900, 900)

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  par(oma=c(2,2,2,2))

  plot(x, y1, xlim=c(0,50), ylim=c(0,10), type="l", xaxt="n", yaxt="n",
     col="red", xlab="x", ylab="axis y1 (red)", main="2 axes (Y)")
  axis(1, 0:5*10)
  axis(2, 0:5*2)

  par(new=T)

  plot(x, y2, xlim=c(0,50), ylim=c(0,1000), type="l", xaxt="n", yaxt="n",
       col="blue", xlab="", ylab="", lwd=2)
  axis(4)
  mtext("axis y2 (blue)", side=4, line=3)
}

#pdf(file="./scatter4_(2-axis-y)).pdf")
png(file="./scatter4_(2-axis-y)).png")
scatter4()
dev.off()


# - - - - - - - - - - - - - - - - - - - -
# scatter5
# 2軸 (X)

scatter5 <- function() {
  x1 <- c(0.1, 0.1, 0.2, 0.4, 0.8, 0.8)
  x2 <- c(420, 300, 120, 130, 50, 20)
  y  <- c(0, 2, 4, 6, 8, 10)

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  par(oma=c(0,1,2,1))

  plot(x1, y, xlim=c(0,1), ylim=c(0,10), type="l", xaxt="n", yaxt="n",
       col="red", xlab="axis x1 (red)", ylab="y", main="2 axes (X)")
  axis(1, 0:5*0.2)
  axis(2, 0:5*2)
  
  par(new=T)

  plot(x2, y, xlim=c(0,500), ylim=c(0,10), type="l", xaxt="n", yaxt="n",
       col="blue", xlab="", ylab="", lwd=2)
  axis(3)
  mtext("axis x2 (blue)", side=3, line=3)
}

#pdf(file="./scatter5_(2-axis-x)).pdf")
png(file="./scatter5_(2-axis-x)).png")
scatter5()
dev.off()




