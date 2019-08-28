demo_deledalle_SR_EPLL_GGMM.m compare the EPLL with the GMM, the LMM and the GGMM. 
It compiles for 5 images with size 512*512. You can change the images for super resolution. 
You can specified the different magnification factor (Q in the code).
The number of components is fixed to 200.
You can choose the size of the patches (tau*tau).
It uses the folders with "EPLL" in their name.
It saves one table (for all Q) in the folder "tab" and the images in the folders "png" and "fig" with a good name.


demo_deledalle_SR_FEPLL.m  compare the EPLL and the FEPLL.
It compiles for 5 images with size 512*512. You can change the images for super resolution. 
You can specified the different magnification factor (Q in the code).
The number of components is fixed to 200.
You can choose the size of the patches (tau*tau).
It uses the folders with "EPLL" and "FEPLL" in their name.
It saves one table (for all Q) in the folder "tab" and the images in the folders "png" and "fig" with a good name.

demo_deledalle_SR_EM_MMSE_EPLL.m compare the MMSE and the EPLL for super resolution.
It compiles 10 times to have stable results and for 5 images with size 512*512. You can change this number ("moyenne" in the code) 
and the images. 
You can specified the different magnification factor (Q in the code).
You can choose the size of the patches (tau*tau).
It uses another EM algorithm then you can specified the different number of components (K in the code).
It uses the folders with "EPLL" in their name and the folders : "Sandeep", "patches", "operators".
It saves one table (for all Q) by part of the image and one image by part of the image and by number Q. 
There are 5 parts: - large: all the image
                   - parthg : top-left part
                   - parthd : top-right part
                   - partbg : bottom-left part
                   - partbd : bottom-right part
In main_sandeep_SR_EM_MMSE_EPLL.m you can choose if you want compute the MMSE using all the components.                  
                   
demo_sandeep_noise_EM_MMSE_EPLL.m compare the MMSE and the EPLL for super resolution.       
It compiles 10 times to have stable results and for 5 images with size 512*512. You can change this number ("moyenne" in the code) 
and the images. 
You can choose the size of the patches (tau*tau).
You can choose the number "sigma" the variance of noise in the image.
It uses another EM algorithm then you can specified the different number of components (K in the code).
It uses the folders with "EPLL" in their name and the folders : "Sandeep", "patches", "operators".


demo_niknejad_SR.m use the local method for image super resolution.
You can change the number of iterations "d", the magnification factor "q" and the size of patches "tau".
You can change the image and the others parameters are well defined for this work.
It saves one table for one q and it saves the reconstructed image.
