/*
This file contains the functions you will modify in order to complete the labs.
These functions allow you to perform DSP on a per-frame or per-sample basis.
When processing on a per-sample basis, there are separate functions for left and right channels.
When processing on a per-frame basis, the left and right channels are interleaved.
The sample rate and frame size can be modified in lab.h.
You can also configure which of the four functions are active in lab.h
*/

#include "main.h"
#include "lab.h"

//

//These functions allow estimation of the number of elapsed clock cycles
extern void tic(void);
extern uint32_t toc(void);

//variables used for the spectrum visualization
extern arm_rfft_fast_instance_f32 fft_inst;
extern float32_t fft_in[FRAME_SIZE/4];
extern float32_t fft_out[FRAME_SIZE/4];
extern float32_t fft_mag[FRAME_SIZE/8];

//declare variables local to this file
uint32_t elapsed_cycles;

/*
This function will be called once before beginning the main program loop.
This is the best place to build a lookup table.
*/
uint32_t tree[512] = {
0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x400,0x0,0x0,0x0,0x1c600,
0x0,0x0,0x0,0x1e100,0x0,0x0,0x0,0x1e000,0x0,0x0,0x0,0x1f080,0x0,0x0,0x0,0x3880,
0x0,0x0,0x0,0x1880,0x0,0x0,0x0,0xc80,0x0,0x0,0x0,0x680,0x0,0x0,0x0,0x680,
0x0,0x0,0x0,0x4680,0x0,0x0,0x0,0x6680,0x0,0x0,0x0,0x680,0x0,0x0,0x0,0xc0e80,
0x0,0x0,0x0,0x40e80,0x0,0x0,0x0,0x20680,0x0,0x0,0x0,0x10680,0x0,0x0,0x0,0x8680,
0x0,0x0,0x0,0x680,0x0,0x0,0x0,0x680,0x0,0x0,0x0,0x680,0x0,0x0,0x0,0x680,
0x0,0x0,0x0,0x600680,0x0,0x0,0x0,0x104680,0x0,0x0,0x0,0x82680,0x0,0x0,0x0,0x41680,
0x0,0x1c00000,0x0,0x20e80,0x0,0x1e00000,0x0,0x10680,0x0,0x3e00000,0x0,0x8680,0x30000000,0x83c00004,0x1,0x4680,
0x20000000,0x8600200c,0x1,0x680,0x40000000,0xc002008,0x0,0x680,0x8f000000,0x18404010,0x0,0x680,0xff000000,0x1c408000,0x0,0x380,
0xff000000,0x9801038c,0x0,0x802180,0x86000000,0x78004399,0x0,0x4010c0,0xc00000,0x30004793,0x0,0x200860,0x8c00000,0x63818fa6,0x21,0x100460,
0x19008000,0x67011c2c,0x71,0x80260,0x22018c00,0x66081838,0xd2,0x40160,0xfff20e00,0x6c0c1f3f,0xf4,0x200e0,0x21b00,0x78021a31,0x68,0x810c60,
0x80009e00,0x7f82d030,0x22,0x38c08c60,0x9c01bc00,0x60807030,0x62,0x3c004260,0x3c626000,0xc0383230,0x6c,0x3c202160,0x3c40e040,0x87f83230,0xc9,0x3c1020e0,
0x7880c040,0xff86260,0x183,0x6082060,0xe100cc80,0x1c39c1ec,0x18e,0x3ff2060,0xfe039900,0xf801818d,0x18f,0x1ff9060,0x80073200,0xf70f0303,0x19f,0x1c860,
0x1e60180,0x11e0603,0x1b0,0x26463,0x20c0180,0x300c06,0x800001fe,0x4363f,0xfffffc00,0x7ffffff,0xe00001c8,0x201f9f,0xfffffc00,0x4fffffff,0xf0000184,0x400fc3,
0x3800180,0xd8030002,0x38000383,0x18000e0,0x11fc0100,0x3fc1f841,0x1c000670,0x1801e60,0xb007c200,0x6200ffe0,0x8ffffe00,0x21ff,0x4003e360,0xc1082707,0xe0000000,0xc01f,
0x23026110,0xe8fc1bff,0xffffffff,0x18007,0x13b3088,0x1c0019ff,0x3fffffff,0x3,0xfed88c,0xc38007,0x0,0x600,0x1fe4c04,0xe0ff80c4,0x1fffffff,0x100603,
0x3f180e00,0x71ff84c4,0xf0000000,0x18010f,0x7e010bc0,0x3be30808,0xe3fff800,0x400bf,0xc41089f0,0x9e0c1e7f,0x87fffe0f,0x38073,0x83f06478,0xe2c46ff,0xc000700,0x43e0,
0x1e00236c,0x1760c1c0,0x180003e0,0x2600,0xfe000338,0x1bffffff,0xb00001ff,0x4001c03,0xe240010,0xf0038007,0xe00000e1,0xf80181f,0xc0000,0xd801c063,0xc0000060,0x40303f,
0xf0080000,0x68f8fe61,0x30,0x207063,0xffffe000,0x3100ff80,0x18,0x1fc060,0x5c011000,0x1a313bb8,0x1c,0xc0060,0x4e00c800,0x1ff919f8,0x1e,0x181060,
0x37f04600,0x7f90cf8,0x1b,0x301060,0xa3fc0600,0x80380639,0x19,0x702060,0xe0000,0xc601f302,0x18,0x3c040e0,0x8b0000,0xf20411c6,0x18,0x7808160,
0x64cf0000,0xf00801e0,0x30,0x3c010660,0x22460000,0xff9f0160,0xf0,0x78220060,0x19000000,0xdcf10130,0x1e0,0x78c23060,0x800000,0xce7f1c1e,0x1e0,0x78c260c0,
0x800000,0xc731fd0f,0x1e0,0x100201c0,0x80000000,0xf310fc89,0x80,0x40300,0xc0000000,0xc9ce18cc,0x0,0x80600,0x78000000,0x66e04002,0x0,0x100e00,
0x78000000,0x74706003,0x0,0x200a00,0x78000000,0x30202000,0x0,0x401a00,0x30000000,0x38000004,0x0,0x803a00,0x0,0x3fb80002,0x0,0x1004a00,
0x80000000,0x37c00001,0x0,0x1008a00,0x80000000,0x30e70001,0x0,0x10a00,0x0,0x307f0000,0x0,0x20a00,0x0,0x703f0000,0x0,0x20a00,
0x0,0x78070000,0x0,0xa00,0x0,0x78000000,0x0,0xa00,0x0,0x70000000,0x0,0x40a00,0x0,0x0,0x0,0x81a00,
0x0,0x0,0x0,0x103a00,0x0,0x0,0x0,0x304a00,0x0,0x0,0x0,0x208a00,0x0,0x0,0x0,0x10a00,
0x0,0x0,0x0,0x20a00,0x0,0x0,0x0,0x40a00,0x0,0x0,0x0,0x80a00,0x0,0x0,0x0,0x100a40,
0x0,0x0,0x0,0x2008c0,0x0,0x0,0x0,0x601800,0x0,0x0,0x0,0xc03a00,0x0,0x0,0x0,0x4a00,
0x0,0x0,0x0,0x38a00,0x0,0x0,0x0,0x40a00,0x0,0x0,0x0,0x80a00,0x0,0x0,0x0,0x101a00,
0x0,0x0,0x0,0x303200,0x0,0x0,0x0,0x207200,0x0,0x0,0x0,0x1e200,0x0,0x0,0x0,0x3c000,
0x0,0x0,0x0,0x3c200,0x0,0x0,0x0,0x18200,0x0,0x0,0x0,0x200,0x0,0x0,0x0,0x400,
0x0,0x0,0x0,0x1800,0x0,0x0,0x0,0x1800,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0
};

//PN-Lab4&Lab5
uint32_t PN[513] = {699842561,3256968010,555344952,788776305,437171664,3911004429,75469212,3644864168,
		1390569626,1204362033,1587358809,1536685593,3152678175,1997803969,3293657319,1583845149,
		2703353362,962108787,2248695597,2400044539,1603822544,1310493618,503032993,1337575995,
		3073722332,3557780031,478222363,3791514872,2550934038,516871296,2716551081,158387105,
		1896591605,3110745373,2315405977,852311565,302096361,3145335264,47169998,2353481449,
		168146218,2037500859,1251152686,610060212,3190760390,711531214,1867498340,3464614586,
		969786542,2130389497,1475916656,507436416,1091312148,2852408006,2879321814,1555542047,
		60126753,1360468804,2212662727,907894519,3816674554,1696231758,3897741930,3985463133,
		321596499,1251798088,3005557977,252328520,337942118,378406773,2195548917,4082706483,
		286989307,2974096497,2876309786,2460171128,3266036827,4145628214,2512158951,657658577,
		361638025,725157733,2823738239,1721826140,1562925357,2039780915,2664879615,3299108037,
		866677883,2461429744,4127809335,3029689888,3266220409,2050089709,90146914,855086774,
		3268318342,4010635103,3997005579,3154176674,4282596397,357686714,800312027,482418353,
		131406372,3585114200,362054503,2296543076,1525147879,268396138,2045472000,622431895,
		1301063857,3250852158,1555582841,2498534732,3761373898,3182362991,1339857638,2938203317,
		2945359315,2276647016,187758684,2791633097,2570308900,3683295812,1401784488,3710887582,
		841507717,3326368056,1535287953,2403898995,882389975,2470652281,56983329,2985314041,
		546704544,3968342228,111170040,3968683510,2341009187,2087601523,2655810372,2871982460,
		1517820798,3982318856,2161642112,1328953352,3961340169,4266862380,880294781,3392031091,
		1937691861,699634979,1326040465,2986011837,986136854,3444443615,757684535,1108096488,
		65021782,4035513754,722227384,4026020733,2213710957,2409063232,806551145,667739319,
		1645256281,963157900,1284905460,4206672463,2565070298,3503972327,1624634835,2944308435,
		1733672405,2823059771,2096215786,2095990822,1377672956,811052513,1288484016,3209641106,
		670452748,3000997174,1830809907,2690953009,3405452804,1839551709,3170519526,99953702,
		2548256214,3736044599,1469965229,103479017,982380024,825733813,1096016498,1658809024,
		2221436805,729218729,3642740430,2587445404,1748992610,693920495,2166181622,1891917511,
		2630426212,3675260159,1007097105,3036662171,1304646013,2962558095,2686581945,2717599235,
		2969237014,2724502118,4223443172,4196010,3861444316,1304209998,190470119,688468173,
		1402413894,2129974943,3229702173,2937832462,2135004860,3552939735,3888271001,2518899842,
		1819624483,3081417797,917328221,3853680539,2345818727,966640558,559538834,3375051693,
		1471291294,3796330218,762749265,2329351662,739687941,2269659948,4052150359,618029221,
		1747944392,2429968216,1383228517,840092478,378229335,267766830,1673142454,81764252,
		1711644798,1866819872,3570441996,403428773,1433542966,4185508718,2522222069,3454751485,
		524577082,4014488791,2235377420,1633634409,3789418413,3944239480,3088089511,608595034,
		502551495,3636470102,101736018,3930252951,1698749450,2913031351,1385744011,1898899074,
		1196877992,3158251323,2728694169,3519396736,3599104121,4224061173,916978457,4293079597,
		2853456748,303352673,2413337228,1094968280,3679803255,1460094230,1772120912,1398637309,
		1323073613,931313605,2575054056,2866528245,2937338728,3896016337,1658885977,3650077233,
		4024604819,539622508,2101237976,1498243935,51368090,1050511101,3929203394,1240886287,
		2511110477,2662936422,3336634906,1777308316,571252556,2990879629,43347082,3384066100,
		3037734140,460587370,1197021134,726072406,327256087,4018681640,2940484712,3080369850,
		4231206788,2416115759,1286202477,2807764987,1551350752,703360325,2278745495,2651616750,
		1292242848,398766384,3859221231,2849061453,480320334,2464782230,1053982513,2607021939,
		3044007887,2696632834,3463568314,3643447571,3717179550,2372793683,4225980546,3365611180,
		3386797967,981329144,3520130824,3794232085,3091868899,1639312007,2732704273,3133533063,
		198859442,3847806325,3371071451,1749055825,2459122055,137724546,2192063362,1388789485,
		3119892100,1755345403,3422443917,2016291190,451668056,1029240796,3834065645,1422352661,
		3929873828,3733948258,612997315,2700387278,3206596107,632498899,3902305233,3718578575,
		3827876747,2096831662,1723794832,1943769079,467576110,3796749486,925363431,2870019301,
		127907530,701053746,2045655074,2821793868,3715080244,3556808025,2348348278,1356395894,
		2789320382,1729040529,3988258790,3286509372,518750967,1607965212,1062049283,3781304161,
		2058228084,3892163328,3418681654,2833326617,1318880999,3517595929,3569393318,2709645330,
		2250050469,3148832407,471089606,420387372,1090083997,1600678608,2935739919,3197054406,
		2502724632,1390985950,1570264711,2136685906,1881079510,359363841,4286366132,1219915790,
		3945245921,3858181139,3771861600,3837944165,1070018834,929218415,3222361314,3670021121,
		928508082,2275914976,1066558256,3857240799,3458325690,2257080952,1130113091,3303094093,
		1487414396,1340722491,3902246071,1255992034,1426971141,1119260678,528400769,1065826232,
		3511306675,2376990380,3516548070,509532287,3559954342,1097250223,627673794,1631197364,
		4190749892,2425988756,589527508,3465246974,597461272,1597847794,2087828415,2965349790,
		3378825825,652411439,517920042,409041438,3659522354,870130956,858290478,1592230472,
		1836059562,2720123374,2428301539,3709838433,4172686684,3019187852,2633571483,294535718,
		1236691621,1938061201,3540354344,3451606525,1088856657,1906031626,3434210578,3256104087,
		1722745914,3394351680,3367079869,2694535511,3177554132,2146749492,1478656877,2582211893,
		1378404641,2011254243,2845754241,1844446975,1852843094,197502078,3416984495,2862892230};



 //LAB5-week1-coefficient
uint32_t data_stream[513] = {0};
float32_t pulse_shaping_coeffs[64] = {
    0.000000f, 0.000163f,-0.000000f,-0.000576f,-0.001624f,-0.003164f,-0.005162f,-0.007522f,
    -0.010081f,-0.012608f,-0.014806f,-0.016326f,-0.016777f,-0.015752f,-0.012850f,-0.007702f,
     0.000000f, 0.010474f, 0.023826f, 0.040029f, 0.058914f, 0.080159f, 0.103297f, 0.127730f,
     0.152748f, 0.177557f, 0.201320f, 0.223195f, 0.242375f, 0.258132f, 0.269854f, 0.277080f,
     0.279521f, 0.277080f, 0.269854f, 0.258132f, 0.242375f, 0.223195f, 0.201320f, 0.177557f,
     0.152748f, 0.127730f, 0.103297f, 0.080159f, 0.058914f, 0.040029f, 0.023826f, 0.010474f,
     0.000000f,-0.007702f,-0.012850f,-0.015752f,-0.016777f,-0.016326f,-0.014806f,-0.012608f,
    -0.010081f,-0.007522f,-0.005162f,-0.003164f,-0.001624f,-0.000576f,-0.000000f, 0.000163f
};

//LAB5-week1-coefficient
arm_fir_interpolate_instance_f32 filter_instance;
uint8_t upsampling_factor = 16;
uint16_t num_taps = 64;
float32_t state[(FRAME_SIZE/4) + 4 - 1] = {0};
float32_t filter_in[FRAME_SIZE/64] = {0};
float32_t filter_out[FRAME_SIZE/4] = {0};
float32_t cos_lut[4] = {1,0,-1,0};
uint32_t i_lut = 0;


/*
void lab_init(void)
{
	  data_stream[0] = 0x967c6ea1;   // header

	    for(int i_word = 0; i_word < 512; i_word++)
	    {
	        //data_stream[i_word+1] = tree[i_word] ^ PN[i_word]; // XOR scrambling
	    	data_stream[i_word+1] = tree[i_word];
	    }
	arm_fir_interpolate_init_f32 (&filter_instance, upsampling_factor, num_taps, pulse_shaping_coeffs, state, FRAME_SIZE/4);

}
*/



void lab_init(void)
{


}

/*
 * 1. DSP Filter Buffer Initialization
    memset(Ux, 0, sizeof(Ux));
    memset(Lx, 0, sizeof(Lx));

    // 2. Costas Loop Variable Initialization
    theta = 0.0f;

    // 3. Symbol Timing & Synchronization Counter Initialization
    k = 0;

    // 4. Data Collection Buffer and Variables Initialization
    i_word = 0;
    i_bit = 0;
    header_matched = 0;
    R = 0;

    memset(data, 0, sizeof(data));
    memset(xcorr, 0, sizeof(xcorr));

    // 5. Baseband Buffer Initialization (Optional, depending on use)
    i_baseband = 0;
    // memset(baseband, 0, sizeof(baseband));
 */


/*
This function will be called each time a complete frame of data is recorded.
Modify this function as needed.
Default behavior:
	1. Deinterleave the left and right channels
	2. Combine the two channels (by addition) into one signal
	3. Save the result to the fft_in buffer which will be used for the display
	4. The original audio buffer is left unchanged (passthrough)
*/
/*

 *

 */

/*
void process_input_buffer(int16_t* input_buffer)
{
	int16_t left_sample;
	int16_t right_sample;
	for (uint32_t i_sample = 0; i_sample<FRAME_SIZE/64; i_sample+=1)
	 {
	     filter_in[i_sample] = ( data_stream[i_word] & (1<<i_bit) ) >> i_bit;
	     filter_in[i_sample] = (filter_in[i_sample]*2) - 1;

	     i_bit += 1;
	     if (i_bit == 32)
	     {
	         i_word += 1;
	         i_bit = 0;
	         if (i_word > 512){i_word = 0;}
	     }
	 }

	arm_fir_interpolate_f32 (&filter_instance, filter_in, filter_out, FRAME_SIZE/64);

		 for (uint32_t i_sample = 0; i_sample < FRAME_SIZE/2; i_sample+=1)
		 {
		      input_buffer[i_sample] = OUTPUT_SCALE_FACTOR*filter_out[i_sample/2]*cos_lut[i_lut];
		      i_lut = (i_lut + 1) % 4;
		      i_sample+=1;
		      input_buffer[i_sample] = 0;
		 }

	return;
}

*/
void process_input_buffer(int16_t* input_buffer)
 {
	int16_t left_sample;
		 	int16_t right_sample;
		 	for (uint32_t i_sample = 0; i_sample < FRAME_SIZE/2; i_sample+=1)
		 	{
		 		left_sample = input_buffer[i_sample];
		 		i_sample +=1;
		 		right_sample = input_buffer[i_sample];
		 		fft_in[i_sample/2] =  (((float32_t) left_sample) + ((float32_t) right_sample))/2;
		 	}
		 	arm_rfft_fast_f32(&fft_inst, fft_in, fft_out, 0);
		 	arm_cmplx_mag_f32(fft_out, fft_mag, FRAME_SIZE/8);
 	return;
 }


/*
 *
void process_input_buffer(int16_t* input_buffer)
 {

 	return;
 }
 *


*/

/*
This function provides access to each individual sample that is incoming on the left channel.
The returned value will be sent to the left channel DAC.
Default behavior:
	1. Copy input to output without modification (passthrough)
	2. Estimate the number of cycles that have elapsed during the function call

*/



/*
 //original process_left
int16_t process_left_sample(int16_t input_sample)
  {

  	tic();
  	int16_t output_sample;
  	output_sample = input_sample;
  	elapsed_cycles = toc();
  	return output_sample;
  }

 */
//LAB5-week2-coefficient
uint32_t i_sample = 0;
float32_t Ux[64] = {0};
float32_t Lx[64] = {0};
float32_t Uy = 0.0f;
float32_t Ly = 0.0f;
float32_t wc = 1.5707963267949f;
float32_t theta = 0.0f;
float32_t mu = 0.01f;
float32_t baseband[FRAME_SIZE/4] = {0};
uint32_t i_baseband = 0;
uint32_t k = 0;
uint32_t data[512] = {0};

//LAB5-week3-coefficeint
int8_t xcorr[32] = {0};
int8_t header[32] = {1,-1,-1,1,-1,1,1,-1,-1,1,1,1,1,1,-1,-1,-1,1,1,-1,1,1,1,-1,1,-1,1,-1,-1,-1,-1,1};
int8_t R = 0;
int8_t header_matched = 0;

//lab5-week1

uint32_t i_word = 0;
uint32_t i_bit = 0;
/*week3
 *
//week2
int16_t process_right_sample(int16_t input_sample)
  {
	        tic();

	 	  	float32_t r = input_sample*INPUT_SCALE_FACTOR;

	 	  		 Ux[0] = r*arm_cos_f32(theta);
	 	  		 Lx[0] = r*arm_sin_f32(theta);

	 	  		 Uy = 0;
	 	  		 Ly = 0;

	 	  		 for (uint32_t i_sample = 0; i_sample < 64; i_sample +=1)
	 	  		 {
	 	  		     Uy += Ux[i_sample]*pulse_shaping_coeffs[i_sample];
	 	  		     Ly += Lx[i_sample]*pulse_shaping_coeffs[i_sample];
	 	  		 }

	 	  		 baseband[i_baseband] = Uy;
	 	  		 i_baseband += 1;
	 	  		 if (i_baseband > FRAME_SIZE/4){i_baseband = 0;}

	 	  		 for (i_sample = 63; i_sample > 0; i_sample -=1)
	 	  		 {
	 	  		     Ux[i_sample] = Ux[i_sample-1];
	 	  		     Lx[i_sample] = Lx[i_sample-1];
	 	  		 }
	 	  		theta += wc - mu*Uy*Ly;
	 	  		if (theta > 6.283185){theta -= 6.283185;}
//
	 	  		if (++k == 16)     // 16 samples per symbol
	 	  		{
	 	  		    k = 0;

	 	  		    if (Uy > 0.0f)
	 	  		        data[i_word] |= (1u << i_bit);

	 	  		    i_bit++;

	 	  		    if (i_bit == 32)
	 	  		    {
	 	  		        i_bit = 0;
	 	  		        i_word++;

	 	  		        if (i_word == 512)
	 	  		        {
	 	  		            display_image(data,128,128);
	 	  		            memset(data,0,sizeof(data));
	 	  		            i_word = 0;
	 	  		        }
	 	  		    }
	 	  		}
  return OUTPUT_SCALE_FACTOR*Uy*0.3;
  }
 *
 */





int16_t process_right_sample(int16_t input_sample)
 {

	 	  	     float32_t r = input_sample*INPUT_SCALE_FACTOR;

	 	  		 Ux[0] = r*arm_cos_f32(theta);
	 	  		 Lx[0] = r*arm_sin_f32(theta);

	 	  		 Uy = 0;
	 	  		 Ly = 0;

	 	  		 for (uint32_t i_sample = 0; i_sample < 64; i_sample +=1)
	 	  		 {
	 	  		     Uy += Ux[i_sample]*pulse_shaping_coeffs[i_sample];
	 	  		     Ly += Lx[i_sample]*pulse_shaping_coeffs[i_sample];
	 	  		 }
	 	  		 //save baseband buffer
	 	  		 baseband[i_baseband] = Uy;
	 	  		 i_baseband += 1;
	 	  		 if (i_baseband > FRAME_SIZE/4){i_baseband = 0;}

	 	  		 for (i_sample = 63; i_sample > 0; i_sample -=1)
	 	  		 {
	 	  		     Ux[i_sample] = Ux[i_sample-1];
	 	  		     Lx[i_sample] = Lx[i_sample-1];
	 	  		 }
	 	  		//costa loops
	 	  		theta += wc - mu*Uy*Ly;
	 	  		if (theta > 6.283185){theta -= 6.283185;}

	 	  		//16 sample down sampling
	 	  		if (++k == 16){
	 	  		        k = 0;
	 	  		        //updating xcorr
	 	  		        if (Uy > 0){xcorr[0] = 1;}
	 	  		        else{xcorr[0] = -1;}
	 	  		     //header cross correlation
	 	  		        R = 0;
	 	  		        for (uint32_t i_hdr = 0; i_hdr < 32; i_hdr +=1){
	 	  		            R += xcorr[i_hdr]*header[i_hdr];}
	 	  		        for (uint32_t i_hdr = 31; i_hdr > 0; i_hdr -=1){
	 	  		            xcorr[i_hdr] = xcorr[i_hdr - 1];}
	 	  		        //maximum value
	 	  		        if (R > 15) {header_matched = 1;}
	 	  		    // Move the data collection code block into an if statement
	 	  		    //so that it only occurs after the header is detected
	 	  		        	if (header_matched) {
	 	  		        		if (Uy > 0.0f)
	 	  		        		data[i_word] |= (1u << i_bit);

	 	  		        	    i_bit++;

	 	  		        	    if (i_bit == 32)
	 	  		        		 {
	 	  		        			 i_bit = 0;
	 	  		        			  i_word++;

	 	  		        		if (i_word == 512)
	 	  		        		{
	 	  		        		display_image(data,128,128);
	 	  		        		memset(data,0,sizeof(data));
	 	  		        		i_word = 0;
	 	  		                   header_matched = 0;
	 	  		                 }
	 	  		        		 }
	 	  		               }
	 	  		            } else {

	 	  		           }

	 return (int16_t)(OUTPUT_SCALE_FACTOR * Uy * 0.3f);
	}


/*week3 circular buffer
 *  float32_t b[31] = {
  0.00346062343384146f,  0.00375345633229333f, -0.00165978419324218f,
  0.01604437543138730f,  0.01160234522751270f, -0.05583228681544590f,
 -0.06588333513764080f,  0.05113597655314760f,  0.09338067914532330f,
  0.00121393293688455f,  0.01031887875865170f,  0.05060955985280240f,
 -0.14532328226207600f, -0.25741269988867300f,  0.09410387502763010f,
  0.38583327674859100f,  0.09410387502763010f, -0.25741269988867300f,
 -0.14532328226207600f,  0.05060955985280240f,  0.01031887875865170f,
  0.00121393293688455f,  0.09338067914532330f,  0.05113597655314760f,
 -0.06588333513764080f, -0.05583228681544590f,  0.01160234522751270f,
  0.01604437543138730f, -0.00165978419324218f,  0.00375345633229333f,
  0.00346062343384146f
};

float32_t x[31] = {0};
float32_t y = 0;
int idx = 0;

int16_t process_left_sample(int16_t input_sample)
{   tic();
   x[idx] = ((float32_t)input_sample) * INPUT_SCALE_FACTOR;
   y = 0;
   int k = idx;
   for (int i = 0; i < 31; i++) {
       y += b[i] * x[k];
       if (--k < 0) k = 30; // wrap around
   }
   if (++idx >= 31) idx = 0; // advance and wrap

	elapsed_cycles = toc();
	return y*OUTPUT_SCALE_FACTOR;
}








 */



/*
This function provides access to each individual sample that is incoming on the left channel.
The returned value will be sent to the right channel DAC.
Default behavior:
	1. Copy input to output without modification (passthrough)
	2. Estimate the number of cycles that have elapsed during the function call
*/
int16_t process_left_sample(int16_t input_sample)
	{
	tic();
			int16_t output_sample;
			output_sample = input_sample;
			elapsed_cycles = toc();
			return output_sample;

	}

/*
This function provides another opportunity to access the frame of data
The default behavior is to leave the buffer unchanged (passthrough)
The buffer you see here will have any changes that occurred to the signal due to:
	1. the process_input_buffer function
	2. the process_left_sample and process_right_sample functions
*/
void process_output_buffer(int16_t* output_buffer)
	{

return;
	}

