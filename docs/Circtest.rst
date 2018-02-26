*******************************************************************************************************************
CircTest module
*******************************************************************************************************************

Prerequisites:

- The `CircTest <https://github.com/dieterich-lab/CircTest>`_ package must be installed


Import of DCC output files into R:
==================================

Using user-generated data
---------------------------

.. code-block:: R

  library(CircTest)

  CircRNACount <- read.delim('CircRNACount',header=T)
  LinearCount <- read.delim('LinearCount',header=T)
  CircCoordinates <- read.delim('CircCoordinates',header=T)

  CircRNACount_filtered <- Circ.filter(circ = CircRNACount,
                                       linear = LinearCount,
                                       Nreplicates = 6,
                                       filter.sample = 6,
                                       filter.count = 5,
                                       percentage = 0.1
                                      )

  CircCoordinates_filtered <- CircCoordinates[rownames(CircRNACount_filtered),]
  LinearCount_filtered <- LinearCount[rownames(CircRNACount_filtered),]

Alternatively, the pre-processed Westholm et al. data from CircTest package may be used:
-----------------------------------------------------------------------------------------

.. code-block:: R

  library(CircTest)

  data(Circ)
  CircRNACount_filtered <- Circ
  data(Coordinates)
  CircCoordinates_filtered <- Coordinates
  data(Linear)
  LinearCount_filtered <- Linear

Test for host-independently regulated circRNAs
====================================================================

Execute the test  
-------------------------------------------------------------------

.. code-block:: R

 test = Circ.test(CircRNACount_filtered,
                  LinearCount_filtered,
                  CircCoordinates_filtered,
                  group=c(rep(1,6),rep(2,6),rep(3,6))
                  )

 # Significant result may be shown in a summary table
 View(test$summary_table)

Visualisation of significantly, host-independently regulated circRNAs
-----------------------------------------------------------------------

.. code-block:: R

 for (i in rownames(test$summary_table))  {
  Circ.ratioplot(CircRNACount_filtered,
                 LinearCount_filtered,
                 CircCoordinates_filtered,
                 plotrow=i,
                 groupindicator1=c(rep('1days',6),rep('4days',6),rep('20days',6)),
                 lab_legend='Ages'
                 )
 }

For further details on the usage of CircTest please refer to the corresponding GitHub project.
