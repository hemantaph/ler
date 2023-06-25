# LeR

`LeR` is a statistical based python package whose core function is to calculate detectable rates of both lensing and unlensed GW events. This calculation very much dependent on the other functionality of the package, which can be subdivided into three parts; 1. Sampling of compact binary source properties, 2. Sampling of lens galaxy characteristics and 3. Solving the lens equation to get image properties of the source. The package as a whole relies on `numpy` array operation, `scipy` interpolation and `multiprocessing` functionality of python to increase speed and functionality without compromising on the ease-of-use. The API of `LeR` is structure such that each functionality mentioned stands on this own right for scientific research but also can be used together as needed. Keys features of `LeR` and its dependencies can be summarized as follows,

- Detectable merger rates: 
    * Calculation not only relies on the properties of simulated events but also on detectability provided by the condition of the GW detectors. For this, `LeR` relies on `gwsnr` for the calculation of optimal signl-to-noise ratio (SNR). Due to prowess of `gwsnr`, rate calulation can be done both for present and future detectors with customizable sensitivities. 
    * Merger rates of both the simulated unlensed and lensed events can be calculated and compared. 
- Sampling GW sources:
    * Distribution of source's red-shift is based on the merger rate density of compact binaries, which can be BBH [cite](https://arxiv.org/abs/2306.03827), BNS [cite](https://arxiv.org/abs/2306.03827), primodial black holes (PBHs) [cite](https://arxiv.org/abs/2306.03827) etc. The code is designed to accomodate easy updates or additions of such distribution by the users in the future. 
    * Sampling of BBH masses is done using `gwcosmo` follwing the powerlaw+peak model. Other related properties are sampled from available priors of `bilby`. Each of them can me manually replaced by the user before feeding in for rate computation.
- Sampling of lens galaxies:
    * Lens distribution of follows [cite](https://arxiv.org/abs/2306.03827). It depends on the sampled source red-shifts and also on the optical depth [cite](https://arxiv.org/abs/2306.03827).
    * `LeR` employs Elliptical Power Law model with external shear (EPL+Shear) model for sampling other features of the galaxy, which is available in the `Lenstronomy` package.
    * Rejection sampling is applied on the above samples on condition that whether event is strongly lensed or not.
- Generation of image properties:
    * Source position is sampled from the caustic in the source plane.
    * Sampled lens' properties and source position is fed in `Lenstronomy` to generate properties of the images.
    * Properties like magnification and time-delay is important as it modifies the source signal strength which in turns changes the SNR and detect ability.
    * `LeR` can handle both super-therhold and sub-threshold events in picking detectable events and rate computation.

`LeR` was written to used by both LIGO scientific collaboration and research students for related works in astrophysics. It is currently use in generating detectable lensing events and GW lensing rates with the available information on current and future detectors. The results will predicts the feasibility of various detectors on the detection of such lensing events. Statistics generated from `LeR` will be use in event validation of the ongoing effort to detected lensed gravitational waves. Lastly, `LeR` was design with upgradability in mind to include additional statistics as required by the related research. 

# Installation

Follow the installation instruction at [ler.readthedoc](https://ler.readthedocs.io/en/latest/installation.html)
