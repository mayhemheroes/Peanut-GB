/**
 * MIT License
 * 
 * Copyright (c) 2015-2019 Lior Halphon
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef apu_h
#define apu_h
#include <stdbool.h>
#include <stdint.h>
#include <stddef.h>

#if __clang__
#define UNROLL _Pragma("unroll")
#elif __GNUC__
#define UNROLL _Pragma("GCC unroll 8")
#else
#define UNROLL
#endif

#define CPU_FREQUENCY 0x400000

/* APU ticks are 2MHz, triggered by an internal APU clock. */

typedef struct
{
    int16_t left;
    int16_t right;
} GB_sample_t;

typedef struct
{
    double left;
    double right;
} GB_double_sample_t;

struct apu_s *apu;

typedef void (*GB_sample_callback_t)(struct apu_s *apu, GB_sample_t *sample);

enum GB_CHANNELS {
    GB_SQUARE_1,
    GB_SQUARE_2,
    GB_WAVE,
    GB_NOISE,
    GB_N_CHANNELS
};

typedef enum {
    GB_HIGHPASS_OFF, // Do not apply any filter, keep DC offset
    GB_HIGHPASS_ACCURATE, // Apply a highpass filter similar to the one used on hardware
    GB_HIGHPASS_REMOVE_DC_OFFSET, // Remove DC Offset without affecting the waveform
    GB_HIGHPASS_MAX
} GB_highpass_mode_t;

enum {
    /* Joypad and Serial */
    GB_IO_JOYP       = 0x00, // Joypad (R/W)
    GB_IO_SB         = 0x01, // Serial transfer data (R/W)
    GB_IO_SC         = 0x02, // Serial Transfer Control (R/W)

    /* Missing */

    /* Timers */
    GB_IO_DIV        = 0x04, // Divider Register (R/W)
    GB_IO_TIMA       = 0x05, // Timer counter (R/W)
    GB_IO_TMA        = 0x06, // Timer Modulo (R/W)
    GB_IO_TAC        = 0x07, // Timer Control (R/W)

    /* Missing */

    GB_IO_IF         = 0x0f, // Interrupt Flag (R/W)

    /* Sound */
    GB_IO_NR10       = 0x10, // Channel 1 Sweep register (R/W)
    GB_IO_NR11       = 0x11, // Channel 1 Sound length/Wave pattern duty (R/W)
    GB_IO_NR12       = 0x12, // Channel 1 Volume Envelope (R/W)
    GB_IO_NR13       = 0x13, // Channel 1 Frequency lo (Write Only)
    GB_IO_NR14       = 0x14, // Channel 1 Frequency hi (R/W)
    /* NR20 does not exist */
    GB_IO_NR21       = 0x16, // Channel 2 Sound Length/Wave Pattern Duty (R/W)
    GB_IO_NR22       = 0x17, // Channel 2 Volume Envelope (R/W)
    GB_IO_NR23       = 0x18, // Channel 2 Frequency lo data (W)
    GB_IO_NR24       = 0x19, // Channel 2 Frequency hi data (R/W)
    GB_IO_NR30       = 0x1a, // Channel 3 Sound on/off (R/W)
    GB_IO_NR31       = 0x1b, // Channel 3 Sound Length
    GB_IO_NR32       = 0x1c, // Channel 3 Select output level (R/W)
    GB_IO_NR33       = 0x1d, // Channel 3 Frequency's lower data (W)
    GB_IO_NR34       = 0x1e, // Channel 3 Frequency's higher data (R/W)
    /* NR40 does not exist */
    GB_IO_NR41       = 0x20, // Channel 4 Sound Length (R/W)
    GB_IO_NR42       = 0x21, // Channel 4 Volume Envelope (R/W)
    GB_IO_NR43       = 0x22, // Channel 4 Polynomial Counter (R/W)
    GB_IO_NR44       = 0x23, // Channel 4 Counter/consecutive, Inital (R/W)
    GB_IO_NR50       = 0x24, // Channel control / ON-OFF / Volume (R/W)
    GB_IO_NR51       = 0x25, // Selection of Sound output terminal (R/W)
    GB_IO_NR52       = 0x26, // Sound on/off

    /* Missing */

    GB_IO_WAV_START  = 0x30, // Wave pattern start
    GB_IO_WAV_END    = 0x3f, // Wave pattern end

    /* Graphics */
    GB_IO_LCDC       = 0x40, // LCD Control (R/W)
    GB_IO_STAT       = 0x41, // LCDC Status (R/W)
    GB_IO_SCY        = 0x42, // Scroll Y (R/W)
    GB_IO_SCX        = 0x43, // Scroll X (R/W)
    GB_IO_LY         = 0x44, // LCDC Y-Coordinate (R)
    GB_IO_LYC        = 0x45, // LY Compare (R/W)
    GB_IO_DMA        = 0x46, // DMA Transfer and Start Address (W)
    GB_IO_BGP        = 0x47, // BG Palette Data (R/W) - Non CGB Mode Only
    GB_IO_OBP0       = 0x48, // Object Palette 0 Data (R/W) - Non CGB Mode Only
    GB_IO_OBP1       = 0x49, // Object Palette 1 Data (R/W) - Non CGB Mode Only
    GB_IO_WY         = 0x4a, // Window Y Position (R/W)
    GB_IO_WX         = 0x4b, // Window X Position minus 7 (R/W)
    // Has some undocumented compatibility flags written at boot.
    // Unfortunately it is not readable or writable after boot has finished, so research of this
    // register is quite limited. The value written to this register, however, can be controlled
    // in some cases.
    GB_IO_DMG_EMULATION = 0x4c,

    /* General CGB features */
    GB_IO_KEY1       = 0x4d, // CGB Mode Only - Prepare Speed Switch

    /* Missing */

    GB_IO_VBK        = 0x4f, // CGB Mode Only - VRAM Bank
    GB_IO_BIOS       = 0x50, // Write to disable the BIOS mapping

    /* CGB DMA */
    GB_IO_HDMA1      = 0x51, // CGB Mode Only - New DMA Source, High
    GB_IO_HDMA2      = 0x52, // CGB Mode Only - New DMA Source, Low
    GB_IO_HDMA3      = 0x53, // CGB Mode Only - New DMA Destination, High
    GB_IO_HDMA4      = 0x54, // CGB Mode Only - New DMA Destination, Low
    GB_IO_HDMA5      = 0x55, // CGB Mode Only - New DMA Length/Mode/Start

    /* IR */
    GB_IO_RP         = 0x56, // CGB Mode Only - Infrared Communications Port

    /* Missing */

    /* CGB Paletts */
    GB_IO_BGPI       = 0x68, // CGB Mode Only - Background Palette Index
    GB_IO_BGPD       = 0x69, // CGB Mode Only - Background Palette Data
    GB_IO_OBPI       = 0x6a, // CGB Mode Only - Sprite Palette Index
    GB_IO_OBPD       = 0x6b, // CGB Mode Only - Sprite Palette Data

    // 1 is written for DMG ROMs on a CGB. Does not appear to have an effect.
    GB_IO_DMG_EMULATION_INDICATION   = 0x6c, // (FEh) Bit 0 (Read/Write)

    /* Missing */

    GB_IO_SVBK       = 0x70, // CGB Mode Only - WRAM Bank
    GB_IO_UNKNOWN2   = 0x72, // (00h) - Bit 0-7 (Read/Write)
    GB_IO_UNKNOWN3   = 0x73, // (00h) - Bit 0-7 (Read/Write)
    GB_IO_UNKNOWN4   = 0x74, // (00h) - Bit 0-7 (Read/Write) - CGB Mode Only
    GB_IO_UNKNOWN5   = 0x75, // (8Fh) - Bit 4-6 (Read/Write)
    GB_IO_PCM_12     = 0x76, // Channels 1 and 2 amplitudes
    GB_IO_PCM_34     = 0x77, // Channels 3 and 4 amplitudes
    GB_IO_UNKNOWN8   = 0x7F, // Unknown, write only
};

/* Speed = 1 / Length (in seconds) */
#define DAC_DECAY_SPEED 20000
#define DAC_ATTACK_SPEED 20000


/* Divides nicely and never overflows with 4 channels and 8 (1-8) volume levels */
#ifdef WIIU
/* Todo: Remove this hack once https://github.com/libretro/RetroArch/issues/6252 is fixed*/
#define MAX_CH_AMP (0xFF0 / 2)
#else
#define MAX_CH_AMP 0xFF0
#endif
#define CH_STEP (MAX_CH_AMP/0xF/8)

struct apu_s
{
	bool global_enable;
	uint8_t apu_cycles;

	uint8_t samples[GB_N_CHANNELS];
	bool is_active[GB_N_CHANNELS];

	uint8_t div_divider; // The DIV register ticks the APU at 512Hz, but is then divided
	// once more to generate 128Hz and 64Hz clocks

	uint8_t lf_div; // The APU runs in 2MHz, but channels 1, 2 and 4 run in 1MHZ so we divide
	// need to divide the signal.

	uint8_t square_sweep_countdown; // In 128Hz
	uint8_t square_sweep_calculate_countdown; // In 2 MHz
	uint16_t new_sweep_sample_legnth;
	uint16_t shadow_sweep_sample_legnth;
	bool sweep_enabled;
	bool sweep_decreasing;

	struct {
		uint16_t pulse_length; // Reloaded from NRX1 (xorred), in 256Hz DIV ticks
		uint8_t current_volume; // Reloaded from NRX2
		uint8_t volume_countdown; // Reloaded from NRX2
		uint8_t current_sample_index; /* For save state compatibility,
						 highest bit is reused (See NR14/NR24's
						 write code)*/

		uint16_t sample_countdown; // in APU ticks (Reloaded from sample_length, xorred $7FF)
		uint16_t sample_length; // From NRX3, NRX4, in APU ticks
		bool length_enabled; // NRX4

	} square_channels[2];

	struct {
		bool enable; // NR30
		uint16_t pulse_length; // Reloaded from NR31 (xorred), in 256Hz DIV ticks
		uint8_t shift; // NR32
		uint16_t sample_length; // NR33, NR34, in APU ticks
		bool length_enabled; // NR34

		uint16_t sample_countdown; // in APU ticks (Reloaded from sample_length, xorred $7FF)
		uint8_t current_sample_index;
		uint8_t current_sample; // Current sample before shifting.

		int8_t wave_form[32];
		bool wave_form_just_read;
	} wave_channel;

	struct {
		uint16_t pulse_length; // Reloaded from NR41 (xorred), in 256Hz DIV ticks
		uint8_t current_volume; // Reloaded from NR42
		uint8_t volume_countdown; // Reloaded from NR42
		uint16_t lfsr;
		bool narrow;

		uint16_t sample_countdown; // in APU ticks (Reloaded from sample_length)
		uint16_t sample_length; // From NR43, in APU ticks
		bool length_enabled; // NR44

		uint8_t alignment; // If (NR43 & 7) != 0, samples are aligned to 512KHz clock instead of
		// 1MHz. This variable keeps track of the alignment.

	} noise_channel;

	bool skip_div_event;
	bool current_lfsr_sample;
	uint8_t apu_registers[0x80];

	struct {
		uint_fast16_t sample_rate;

		double sample_cycles; // In 4 MHz units
		double cycles_per_sample;

		// Samples are NOT normalized to MAX_CH_AMP * 4 at this stage!
		unsigned cycles_since_render;
		unsigned last_update[GB_N_CHANNELS];
		GB_sample_t current_sample[GB_N_CHANNELS];
		GB_sample_t summed_samples[GB_N_CHANNELS];

		double dac_discharge[GB_N_CHANNELS];
		GB_highpass_mode_t highpass_mode;
		double highpass_rate;
		GB_double_sample_t highpass_diff;

		GB_sample_callback_t sample_callback;
	} apu_output;

	void *priv;
};

#if 0
void GB_set_highpass_filter_mode(struct apu_s *apu, GB_highpass_mode_t mode);
#endif
void GB_apu_div_event(struct apu_s *apu);
void GB_set_sample_rate(struct apu_s *apu, const uint_fast16_t sample_rate);
void GB_apu_set_sample_callback(struct apu_s *apu, GB_sample_callback_t callback);
bool GB_apu_is_DAC_enabled(struct apu_s *apu, unsigned index);
void GB_apu_write(struct apu_s *apu, uint8_t reg, uint8_t value);
uint8_t GB_apu_read(struct apu_s *apu, uint8_t reg);
void GB_apu_init(struct apu_s *apu);
void GB_apu_run(struct apu_s *apu);

#endif /* apu_h */
