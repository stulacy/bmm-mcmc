#data("K2_N100_P5")
#foo <- gibbs_full(K2_N100_P5, K=2, nsamples=1000, burnin=100)
#plot_gibbs(foo)
#
## Plot thetas
## Original
#plot_gibbs(foo, cluster_threshold = 0.15)
#
## Relabelled
#bar <- foo
#bar$theta <- bar$theta_relabelled
#bar$z <- bar$z_relabelled
#plot_gibbs(bar, cluster_threshold = 0.15)
#
#new <- gibbs_dp(df_3, nsamples=10000, burnin=1000, burnrelabel=500, maxK=20)
#
## Plot thetas
## Original
#plot_gibbs(new, cluster_threshold = 0.15)
#
## Relabelled
#new_relab <- new
#new_relab$theta <- new_relab$theta_relabelled
#new_relab$z <- new_relab$z_relabelled
#plot_gibbs(new_relab, cluster_threshold = 0.15)
#
#
#res_full <- gibbs_full(df_3, nsamples=10000, burnin=1000, K = 3)
#plot_gibbs(res_full)
#
#res_collapsed <- gibbs_collapsed(df_3, nsamples=10000, burnin=1000, K = 3)
#plot_gibbs(res_collapsed)
