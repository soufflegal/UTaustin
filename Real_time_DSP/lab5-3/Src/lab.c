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
void lab_init(int16_t* output_buffer)
{
	return;
}

/*
This function will be called each time a complete frame of data is recorded.
Modify this function as needed.
Default behavior:
	1. Deinterleave the left and right channels
	2. Combine the two channels (by addition) into one signal
	3. Save the result to the fft_in buffer which will be used for the display
	4. The original audio buffer is left unchanged (passthrough)
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
This function provides access to each individual sample that is incoming on the left channel.
The returned value will be sent to the left channel DAC.
Default behavior:
	1. Copy input to output without modification (passthrough)
	2. Estimate the number of cycles that have elapsed during the function call
*/


float32_t Ux[64] = {0};
float32_t Lx[64] = {0};
float32_t Uy = 0.0f;
float32_t Ly = 0.0f;
float32_t wc = 1.5707963267949f;
float32_t theta = 0.0f;
float32_t mu = 0.01f;
float32_t baseband[FRAME_SIZE/4] = {0};
uint32_t i_baseband = 0;
//uint32_t k = 0;
//uint32_t data[512] = {0};



float32_t B1[3] = {
		 0.012699218037574,
		 0,
		 -0.012699218037574
};
float32_t A1[3] = {
		 1,         -1.93670162948468 ,        0.974601563924852
};

// 상태: x[n-1], x[n-2], y[n-1], y[n-2]
float32_t X1[3] = {0};
float32_t Y1[3] = {0};


// ==== BPF2: Elliptic 2nd-order bandpass ====
// scale factor: 0.0250809879172544
// SOS: [1  0  -1   1  -1.80156978709300   0.949838024165491]

float32_t B2[3] = {
		0.0250809879172544 ,                        0,      -0.0250809879172544
};
float32_t A2[3] = {
		1,           -1.801569787093 ,        0.949838024165491
};

float32_t X2[3] = {0};
float32_t Y2[3] = {0};












//LAB5-week1-coefficient
//uint32_t data_stream[513] = {0};
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

int16_t process_left_sample(int16_t input_sample)
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

		 	  		 for (uint32_t i_sample = 63; i_sample > 0; i_sample -=1)
		 	  		 {
		 	  		     Ux[i_sample] = Ux[i_sample-1];
		 	  		     Lx[i_sample] = Lx[i_sample-1];
		 	  		 }
		 	  		//costa loops
		 	  		theta += wc - mu*Uy*Ly;
		 	  		if (theta > 6.283185){theta -= 6.283185;}
		 	  		// BPF1: 2nd-order IIR (Direct Form I)

		 	  		X1[0]= Uy;
		 	  		Y1[0]=B1[0]*X1[0];
		 	  		for(uint32_t i=1;i<3;i++){

		 	  			Y1[0]+=B1[i]*X1[i] - A1[i]*Y1[i];
		 	  			};
		 	  		for (uint32_t i=2;i>0;i--){
		 	  			X1[i]=X1[i-1];
		 	  			Y1[i]=Y1[i-1];
		 	  		}
		 	  			//Squaring
		 	  			X2[0]=Y1[0]*Y1[0];
		 	  			// second BPF
		 	  			Y2[0]=B2[0]*X2[0];
		 	  			for(uint32_t i=1;i<3;i++){
		 	  			Y2[0]+=B2[i]*X2[i] - A2[i]*Y2[i];
		 	  			};
		 	  			for (uint32_t i=2;i>0;i--){
		 	  			X2[i]=X2[i-1];
		 	  			Y2[i]=Y2[i-1];
		 	  			}



 return Y2[0]*OUTPUT_SCALE_FACTOR*2;
}






/*
This function provides access to each individual sample that is incoming on the left channel.
The returned value will be sent to the right channel DAC.
Default behavior:
	1. Copy input to output without modification (passthrough)
	2. Estimate the number of cycles that have elapsed during the function call
*/
int16_t process_right_sample(int16_t input_sample)
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
