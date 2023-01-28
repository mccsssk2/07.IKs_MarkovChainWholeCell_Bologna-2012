# Repository.  

Codes to estimate parameters of the IKs ion channel under control and PUFA conditions.

# Background.  

This work was done in Bologna, Italy in consultation with Dr. Stefano Severi.
Due to the nature of the workplace, a Matlab program was deemed suitable.
During this work, I co-supervised a student, Dr. Tomas Stary, who progressed
to do a PhD with Vadim Biktashev in Exeter. He is now a data scientist.  

This work is related to the article:  
https://academic.oup.com/cardiovascres/article/105/2/223/546340?login=true  

# Dependencies.

This program uses Matlab 2012+.

# Install.

No install required.

# Sources and data description.

The sources are the matlab files, data is also provided for example.
The matlab files consist of a driver, a RHS of the IKs ion channel Markov chain, and 
the steepest gradient parameter estimator function.

# Use.  

Run the driver.m programs. The program can be (should be) compiled using the matlab compiler
"mcc -mv -R -nodisplay driver.m" to run on clusters/non-interactively. The markov chain 
RHS can be replaced by any other corresponding to data, e.g. INa, ICaL, ...

# Maintainer.  

This code is provided as is.

# Acknowledements.

This project was generously funded by CINECA (Italy), Uni. Bologna (Italy). 

# Licence.  

Code is unsharable as of now.

## Dr. Sanjay R. Kharche, Ph.D. (App. Maths).  
January 23rd, 2023.  
Email: Sanjay.Kharche@lhsc.on.ca , skharche@uwo.ca , srk.london@yahoo.com .  
Phone: +15198780685 (mobile).  
Website: https://kcru.lawsonresearch.ca/research/srk/index.html  

