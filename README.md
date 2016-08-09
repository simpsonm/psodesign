# psomcmc

Some todos:

* Metric: number of iterations until convergence - mean/median & SE or CI
* Check PSO literature for other adaptive approaches to $\omega$ (inertia)

Some notes from Scott that need to be addressed:
* spell check
* pg 6 - are the algorithms sensitive to c?
* does the swarm size impact which target rates are good?
* be careful of switching between past and present tense
* general paren rules: [{()}]
* explain nbhd topologies earlier - before referencing them
* remove the small c speculation - or mention it as a foil to set up "simulations prove this wrong"
* figure out if "tall" data is worth mentioning. If not, there's a mention in the main body of the text
* is it ok to start a sentence with mathematical notation?
* Have 2015 ACS numbers been released? Not ACS, but other pop numbers
* Pg 14 "any set of spatial basis functions" - not true, has to be Areal.
  ****************
  Me: I meant areal, maybe "any set of areal spatial basis functions"?
* pg 16: 'allowing the diagonal elements of \bm{L} to be negative could be problematic'
  ****************
  Me: only problem is people misunderstanding what's going on. 
* pg 16 bottom 'where are the results'
  uh... later
* pg 17 what are the age & educ categories?
  ****************
  can't find it! worthing digging hard for?
* pg 17 is remove the age & education linear effects justifiable?
  ****************
  me: if they were fixed effects they would be unidentified, so kind of redundant.
  * redundant without a sum constraint, or dropping the first of each set of categories, etc.
* pg 20 & 21 - what are these plots telling us?
  ****************
  there's a paragraph describing this.
* pg 22 what does 'minimal tweaking' mean?
  ****************
  me: change adapt_delta and size of search tree - I don't want to get into it, honestly
* remove full details of alternative MCMC algorithms - just cite and say this is how we adapt

An old outline (pre-JSM):
(Note: put all sims except ACC and MCMC in appendix)
* MCMC in complex models, e.g. "tall" data
* normal approximation & IMH/IMHwG - GLM framework and necessary conditions
* problem: finding posterior mode (and conditional mode in IMHwG)
* PSO - PSO, BBPSO(xp)(-MC), neighborhoods
* AT-PSO and AT-BBPSO
* General PSO conclusions
* Application to MCMC for founty population estimates
  ACC sims and MCMC sims
* Discussion - bad for unemployment rates
