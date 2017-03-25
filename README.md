# CV-Playground
My computer vision playground, specifically for __multi-view geometry and reconstruction__ with tested basic Matlab implementation. List here only for reference. All this basic algorithms with better implementation could be found in CV toolboxes like [OpenCV](http://opencv.org/).

### Function documentation
- _calcFSevenPoints.m_: Get the fundamental matrix F from seven correspondences in two given images, using seven points algorithm. [Seven Point Algorithm](http://www.cs.unc.edu/~marc/tutorial/node55.html)

- _calcHFourPoints.m_: Get the projection matrix H from four correspondences in two given images, using four points algorithm. [Four Points Algorithm](http://www.math.kth.se/math/forskningsrapporter/philip.pdf)

- _calcMinorEigenImageAndCorner.m_: Corner detection based gradient image in  five points central difference operator. Forstner corner point
operator is used to extract the feature point.

- _calcSampsonError.m_: Calculate the Sampson error for the re-projection from estimated 3D scene points to 2D image points using H (image projection matrix). See this page [Sampson Error](http://www.ce.unipr.it/people/medici/geometry/node71.html).

- _calcSampsonCorrection.m_: Sampson correction for 3D scene points.

- _ePNP.m_: Implementation of efficient PnP algorithm (A nearly O(n) implementation of PnP algorithm). See this page for details (ePNP)[https://pdfs.semanticscholar.org/7655/247308e345bc6e6b9365461f04b488b56429.pdf].

- _triangulation.m_: Scene points triangulation.

- _umeyama.m_: Implementation of Umeyama's method. See this page for details. (Umeyama's Method)[http://graphics.stanford.edu/courses/cs164-09-spring/Handouts/paper_Umeyama.pdf]
