# Evolution_paper_2024_code
Matlab code used for the computations published in 

Barna Páll-Gergely, András Á. Sipos, Mathias Harzhauser, Aydin Örstan, Viola Winkler, Thomas A. Neubauer (2024): "Many roads to success: Alternative routes to building an economic shell in land snails", Evolution (accepted)

The code uses the following Matlab toolboxes:
(i) geom3D (including mesh3D)
https://nl.mathworks.com/matlabcentral/fileexchange/24484-geom3d?requestedDomain=
(ii) chebfun
http://www.chebfun.org/

The code needs the CT-scanned data of the gastropod shells published at https://doi.org/10.57756/g7jjq4.

The analysis works as follows:
(1) The CT scanned data should be placed in a 'CT' folder placed beside the source code published here.

(2) The source code should include the data.txt file that contains the tilting and translational corrections.

(3) Preprocess.m analyzes the contents of the CT folder, makes the needed corrections, and saves the data in the 'surface' folder and the simplified midsurface to the 'midsurface' folder.

(4) Analysis.m produces the vertical slices and computes all the data needed for the cross-sectional analysis. 

(5) Because in some cases the automatic identification of the cross-section fails, the SectionByHand.m can be used to manually correct the cross-sections. The graphical interface works as follows: 
- points can be placed by the left mouse button
- from point 1 until the key SPACE is pushed the points belong to the outer surface
- points placed after the SPACE is pushed belong to the overlapped surface
- placing points is finished by pushing the ENTER
The corrected cross-sections are stored in the folder SecByHandResults. Note that SecByHandResults0 contains the cross-sections used to produce the figures and tables of the paper. This data is shared upon request (siposa@eik.bme.hu)

(6) AR.m computes the aspect ratio of the shell (later transformed to the epsilon angle in the paper)

(7) IP2D.m writes the 2D isoperimetric ratio into an xls file (for all sections stored in the files of the folder SecByHandResults).

(8) IP3D.m computes the 3D isoperimetric ratio for the selected or all shells. E.g. IP3D([1,2,3]) computes the IP3D for the first three scans, IP3D(-1) for all available STL files.