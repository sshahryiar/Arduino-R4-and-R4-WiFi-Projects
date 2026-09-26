#include <Arduino.h>
#include "Arduino_LED_Matrix.h"
#include <string.h>  // for memset

ArduinoLEDMatrix matrix;

// ===== CONFIGURATION =====
#define N                       64                       // FFT size (must be power of 2)
#define SAMPLE_RATE             8000                     // Sample Rate in Hz
#define MIC_PIN                 A0
#define PERIOD_US               (1000000UL / SAMPLE_RATE)

// Display scaling options
#define USE_AUTO_GAIN           false                   // Automatic level adjustment
#define USE_LOG_SCALE           false                   // Logarithmic frequency binning
#define MANUAL_SCALE_FACTOR     40.0                    // Used when AUTO_GAIN is false
#define MIN_DISPLAY_HEIGHT      0                       // Noise floor threshold

// Pre-calculated constants
#define LOG2_N                  6                       // log2(64) = 6
#define N_HALF                  (N >> 1)                // N/2 = 32
#define BINS_PER_COLUMN         ((N_HALF - 1) / 12)     // Pre-calculate bins per column

// ===== GLOBAL BUFFERS =====
float real[N];
float imag[N];
float mag[N_HALF];
uint8_t frame[8][12];


#if USE_LOG_SCALE
static const uint8_t bin_edges[13] PROGMEM = {1, 2, 3, 5, 7, 10, 14, 19, 25, 33, 43, 56, N_HALF};
#endif


// ===== Bit Reversal =====
static inline uint16_t reverse_bits(uint16_t x, const uint8_t bits) 
{
    uint16_t y = 0;
    for(uint8_t i = 0; i < bits; i++) 
    {
        y = (y << 1) | (x & 1);
        x >>= 1;
    }
    return y;
}


// ===== Hann Window Function =====
void apply_Hann_window(float *data) 
{
    const float factor = ((2.0 * PI) / (N - 1)); 
    
    for(uint8_t i = 0; i < N; i++) 
    {
        data[i] *= (0.5 * (1.0 - cos(factor * i)));
    }
}


// ===== Cooley-Tukey FFT =====
void fft(float *real, float *imag) 
{
    for(uint8_t i = 0; i < N; i++) 
    {
        uint8_t j = reverse_bits(i, LOG2_N);

        if(j > i) 
        {
            // Swap real and imaginary parts
            float temp = real[i];
            real[i] = real[j];
            real[j] = temp;

            temp = imag[i];
            imag[i] = imag[j];
            imag[j] = temp;
        }
    }

    for(uint8_t s = 1; s <= LOG2_N; s++) 
    {
        const uint8_t m = (1 << s);
        const float angle = (-2.0 * PI) / m;
        const float wm_real = cos(angle);
        const float wm_imag = sin(angle);

        for(uint8_t k = 0; k < N; k += m) 
        {
            float w_real = 1.0;
            float w_imag = 0.0;

            const uint8_t half_m = (m >> 1);
            for(uint8_t j = 0; j < half_m; j++) 
            {
                const uint8_t u = k + j;
                const uint8_t t = u + half_m;

                const float tr = (w_real * real[t]) - (w_imag * imag[t]);
                const float ti = (w_real * imag[t]) + (w_imag * real[t]);

                real[t] = real[u] - tr;
                imag[t] = imag[u] - ti;
                real[u] += tr;
                imag[u] += ti;

                const float tmp = w_real;
                w_real = (tmp * wm_real) - (w_imag * wm_imag);
                w_imag = (tmp * wm_imag) + (w_imag * wm_real);
            }
        }
    }
}


// ===== ADC sampling =====
void sample_ADC(void) 
{
    float dc_sum = 0.0;

    for(uint8_t i = 0; i < N; i++) 
    {
        const uint32_t t0 = micros();
        
        const uint16_t raw_adc = analogRead(MIC_PIN);
        real[i] = (float)raw_adc;
        dc_sum += raw_adc;
        imag[i] = 0.0;

        while((micros() - t0) < PERIOD_US);
    }

    const float dc_offset = (dc_sum / N);
    for(uint8_t i = 0; i < N; i++) 
    {
        real[i] -= dc_offset;
    }
}


// ===== Compute magnitude spectrum =====
void compute_magnitude_spectrum(void) 
{
    for(uint8_t i = 0; i < N_HALF; i++) 
    {
        mag[i] = sqrt((real[i] * real[i]) + (imag[i] * imag[i]));
    }
}


// ===== Auto gain computation =====
static inline float calculate_scale(void) 
{
    float scale = 0.0;
    
    if(USE_AUTO_GAIN) 
    {
        float max_mag = 0.0;
        
        for(uint8_t i = 1; i < N_HALF; i++) 
        {
            if(mag[i] > max_mag) 
            {
                max_mag = mag[i];
            }
        }
        
        if(max_mag > 0) 
        {
            scale = (8.0 / max_mag);
        } 
        else 
        {
            scale = (1.0 / MANUAL_SCALE_FACTOR);
        }
    } 
    else 
    {
        scale = (1.0 / MANUAL_SCALE_FACTOR);
    }
    
    return scale;
}


// ===== Display: Linear frequency binning =====
void draw_linear_spectrum(void) 
{
    memset(frame, 0, sizeof(frame));

    const float scale = calculate_scale();

    for(uint8_t col = 0; col < 12; col++) 
    {
        float sum = 0.0;  // Reset for each column

        const uint8_t bin_start = (col * BINS_PER_COLUMN) + 1;
        const uint8_t bin_end = bin_start + BINS_PER_COLUMN;
        
        for(uint8_t bin = bin_start; bin < bin_end && bin < N_HALF; bin++) 
        {
            sum += mag[bin];
        }

        const float avg_level = sum / BINS_PER_COLUMN;
        uint8_t height = constrain((int)(avg_level * scale), 0, 8);

        #if MIN_DISPLAY_HEIGHT > 0
        if(height < MIN_DISPLAY_HEIGHT) 
        {
            height = 0;
        }
        #endif

        for(uint8_t row = 0; row < height; row++) 
        {
            frame[7 - row][col] = 1;
        }
    }

    matrix.renderBitmap(frame, 8, 12);
}


// ===== Display : Logarithmic frequency binning =====
#if USE_LOG_SCALE
void draw_logarithmic_spectrum(void) 
{
    memset(frame, 0, sizeof(frame));

    const float scale = calculate_scale();

    for(uint8_t col = 0; col < 12; col++) 
    {
        float sum = 0.0;
        
        const uint8_t edge_start = pgm_read_byte(& bin_edges[col]);
        const uint8_t edge_end = pgm_read_byte(& bin_edges[col + 1]);
        const uint8_t count = edge_end - edge_start;

        for(uint8_t i = edge_start; i < edge_end; i++) 
        {
            sum += mag[i];
        }

        const float avg_level = (sum / count);
        uint8_t height = constrain((int)(avg_level * scale), 0, 8);

        #if MIN_DISPLAY_HEIGHT > 0
        if(height < MIN_DISPLAY_HEIGHT) 
        {
            height = 0;
        }
        #endif

        for(uint8_t row = 0; row < height; row++) 
        {
            frame[7 - row][col] = 1;
        }
    }

    matrix.renderBitmap(frame, 8, 12);
}
#endif


void setup(void) 
{
    #if (N & (N - 1)) != 0
        #error "N must be a power of 2"
    #endif

    analogReadResolution(14);
    matrix.begin();
}


void loop(void) 
{
    sample_ADC();
    apply_Hann_window(real);
    fft(real, imag);
    compute_magnitude_spectrum();

    #if USE_LOG_SCALE
        draw_logarithmic_spectrum();
    #else
        draw_linear_spectrum();
    #endif
}