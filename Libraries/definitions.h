//#define baseDataFolder "D:\\Data\\Mahmoud Files\\Research\\thesis\\data\\"
//#define baseDataFolder "E:\\mahmoud\\COLLEGE\\Research\\Thesis\\Data"
#define logFileFolder ".\\"
#define logFileName "HyperspectralAnalystQT.log"


//FOR CUDA code
#define CUDA_N (64)
#define MAX_REGION_ADJ (512)
#define HARD_CODED_N_BANDS 220
#define HARD_CODED_1_DIM_CUDA_KERNELS_BLOCK_WIDTH 128

#ifdef AMP_PRECISE_FLOAT
#define amp_math Concurrency::precise_math
#else
#define amp_math Concurrency::fast_math
#endif

#ifdef AMP_PRECISE_FLOAT
#define amp_float float
#else
#define amp_float float
#endif
