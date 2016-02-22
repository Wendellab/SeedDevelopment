# Plot multiple scale-free P(k) lines in one Plot

multiScaleFreePlot<-function(degree, title = "")
{
    n=dim(degree)[2]
    color=brewer.pal(n,"Set1")
    k= degree[,1]
    discretized.k = cut(k, 20)
    dk = tapply(k, discretized.k, mean)
    p.dk = as.vector(tapply(k, discretized.k, length)/length(k))
    plot(log(dk, 10), log(p.dk, 10), col=color[1], main=title)# power law, linear
    ll<-lm(log(p.dk, 10)~log(dk, 10))
    abline(ll, col=color[1])
    legend(min(log(dk, 10)),min(log(p.dk, 10))+0.8, legend=names(degree), pch = 1,  col=color)
    for(j in 2:n)
    {
       k= degree[,j]
       discretized.k = cut(k, 20)
       dk = tapply(k, discretized.k, mean)
       p.dk = as.vector(tapply(k, discretized.k, length)/length(k))
       points(log(dk, 10), log(p.dk, 10), col=color[j])    # power law, # linear
       ll<-lm(log(p.dk, 10)~log(dk, 10))
       abline(ll, col=color[j])
    }
}