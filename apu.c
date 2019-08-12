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

#include <stdint.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include "apu.h"

#define likely(x)   __builtin_expect((x), 1)
#define unlikely(x) __builtin_expect((x), 0)

static const uint8_t duties[] = {
    0, 0, 0, 0, 0, 0, 0, 1,
    1, 0, 0, 0, 0, 0, 0, 1,
    1, 0, 0, 0, 0, 1, 1, 1,
    0, 1, 1, 1, 1, 1, 1, 0,
};

static void refresh_channel(struct apu_s *apu, unsigned index,
		unsigned cycles_offset)
{
    unsigned multiplier = apu->apu_output.cycles_since_render + cycles_offset -
	    apu->apu_output.last_update[index];
    apu->apu_output.summed_samples[index].left +=
	    apu->apu_output.current_sample[index].left * multiplier;
    apu->apu_output.summed_samples[index].right +=
	    apu->apu_output.current_sample[index].right * multiplier;
    apu->apu_output.last_update[index] = apu->apu_output.cycles_since_render +
	    cycles_offset;
}

bool GB_apu_is_DAC_enabled(struct apu_s *apu, unsigned index)
{
#if 0
    if (apu->model >= GB_MODEL_AGB) {
        /* On the AGB, mixing is done digitally, so there are no per-channel
           DACs. Instead, all channels are summed digital regardless of
           whatever the DAC state would be on a CGB or earlier model. */
        return true;
    }
#endif
    
    switch (index) {
        case GB_SQUARE_1:
            return apu->apu_registers[GB_IO_NR12] & 0xF8;

        case GB_SQUARE_2:
            return apu->apu_registers[GB_IO_NR22] & 0xF8;

        case GB_WAVE:
            return apu->wave_channel.enable;

        case GB_NOISE:
            return apu->apu_registers[GB_IO_NR42] & 0xF8;
    }

    return false;
}

static void update_sample(struct apu_s *apu, unsigned index, int8_t value,
		unsigned cycles_offset)
{
#if 0
    if (apu->model >= GB_MODEL_AGB) {
	/* On the AGB, because no analog mixing is done, the behavior of NR51 is
	 * a bit different.
	   A channel that is not connected to a terminal is idenitcal to a
	   connected channel playing PCM sample 0. */
        apu->samples[index] = value;
        
        if (apu->apu_output.sample_rate) {
            unsigned right_volume = (apu->apu_registers[GB_IO_NR50] & 7) + 1;
            unsigned left_volume = ((apu->apu_registers[GB_IO_NR50] >> 4) & 7) + 1;
            
            if (index == GB_WAVE) {
                /* For some reason, channel 3 is inverted on the AGB */
                value ^= 0xF;
            }
            
            GB_sample_t output;
            if (apu->apu_registers[GB_IO_NR51] & (1 << index)) {
                output.right = (0xf - value * 2) * right_volume;
            }
            else {
                output.right = 0xf * right_volume;
            }
            
            if (apu->apu_registers[GB_IO_NR51] & (0x10 << index)) {
                output.left = (0xf - value * 2) * left_volume;
            }
            else {
                output.left = 0xf * left_volume;
            }
            
            if (*(uint32_t *)&(apu->apu_output.current_sample[index]) != *(uint32_t *)&output) {
                refresh_channel(apu, index, cycles_offset);
                apu->apu_output.current_sample[index] = output;
            }
        }
        
        return;
    }
#endif
    
    if (!GB_apu_is_DAC_enabled(apu, index)) {
	    value = apu->samples[index];
    }
    else {
	    apu->samples[index] = value;
    }

    if (apu->apu_output.sample_rate) {
        unsigned right_volume = 0;
        if (apu->apu_registers[GB_IO_NR51] & (1 << index)) {
            right_volume = (apu->apu_registers[GB_IO_NR50] & 7) + 1;
        }
        unsigned left_volume = 0;
        if (apu->apu_registers[GB_IO_NR51] & (0x10 << index)) {
            left_volume = ((apu->apu_registers[GB_IO_NR50] >> 4) & 7) + 1;
        }
        GB_sample_t output = {(0xf - value * 2) * left_volume, (0xf - value * 2) * right_volume};
        if (*(uint32_t *)&(apu->apu_output.current_sample[index]) != *(uint32_t *)&output) {
            refresh_channel(apu, index, cycles_offset);
            apu->apu_output.current_sample[index] = output;
        }
    }
}

static double smooth(double x)
{
    return 3*x*x - 2*x*x*x;
}

static void render(struct apu_s *apu)
{
    GB_sample_t output = {0,0};

    UNROLL
    for (unsigned i = 0; i < GB_N_CHANNELS; i++) {
        double multiplier = CH_STEP;
        
	if (!GB_apu_is_DAC_enabled(apu, i)) {
		apu->apu_output.dac_discharge[i] -= ((double) DAC_DECAY_SPEED) / apu->apu_output.sample_rate;
		if (apu->apu_output.dac_discharge[i] < 0) {
			multiplier = 0;
			apu->apu_output.dac_discharge[i] = 0;
		}
		else {
			multiplier *= smooth(apu->apu_output.dac_discharge[i]);
		}
	}
	else {
		apu->apu_output.dac_discharge[i] += ((double) DAC_ATTACK_SPEED) / apu->apu_output.sample_rate;
		if (apu->apu_output.dac_discharge[i] > 1) {
			apu->apu_output.dac_discharge[i] = 1;
		}
		else {
			multiplier *= smooth(apu->apu_output.dac_discharge[i]);
		}
	}

	if (likely(apu->apu_output.last_update[i] == 0)) {
		output.left += apu->apu_output.current_sample[i].left * multiplier;
		output.right += apu->apu_output.current_sample[i].right * multiplier;
	}
	else {
		refresh_channel(apu, i, 0);
		output.left += (signed long) apu->apu_output.summed_samples[i].left * multiplier
			/ apu->apu_output.cycles_since_render;
		output.right += (signed long) apu->apu_output.summed_samples[i].right * multiplier
			/ apu->apu_output.cycles_since_render;
		apu->apu_output.summed_samples[i] = (GB_sample_t){0, 0};
	}
	apu->apu_output.last_update[i] = 0;
    }
    apu->apu_output.cycles_since_render = 0;

    GB_sample_t filtered_output = apu->apu_output.highpass_mode?
	    (GB_sample_t) {output.left - apu->apu_output.highpass_diff.left,
		    output.right - apu->apu_output.highpass_diff.right} :
			    output;

    switch (apu->apu_output.highpass_mode) {
        case GB_HIGHPASS_OFF:
            apu->apu_output.highpass_diff = (GB_double_sample_t) {0, 0};
            break;
        case GB_HIGHPASS_ACCURATE:
            apu->apu_output.highpass_diff = (GB_double_sample_t)
                {output.left - filtered_output.left * apu->apu_output.highpass_rate,
                    output.right - filtered_output.right * apu->apu_output.highpass_rate};
            break;
	default:
        case GB_HIGHPASS_REMOVE_DC_OFFSET: {
            unsigned mask = apu->apu_registers[GB_IO_NR51];
            unsigned left_volume = 0;
            unsigned right_volume = 0;
            UNROLL
            for (unsigned i = GB_N_CHANNELS; i--;) {
                if (apu->is_active[i]) {
                    if (mask & 1) {
                        left_volume += (apu->apu_registers[GB_IO_NR50] & 7) * CH_STEP * 0xF;
                    }
                    if (mask & 0x10) {
                        right_volume += ((apu->apu_registers[GB_IO_NR50] >> 4) & 7) * CH_STEP * 0xF;
                    }
                }
                else {
                    left_volume += apu->apu_output.current_sample[i].left * CH_STEP;
                    right_volume += apu->apu_output.current_sample[i].right * CH_STEP;
                }
                mask >>= 1;
            }
            apu->apu_output.highpass_diff = (GB_double_sample_t)
            {left_volume * (1 - apu->apu_output.highpass_rate) + apu->apu_output.highpass_diff.left * apu->apu_output.highpass_rate,
                right_volume * (1 - apu->apu_output.highpass_rate) + apu->apu_output.highpass_diff.right * apu->apu_output.highpass_rate};

        case GB_HIGHPASS_MAX:;
        }

    }
    
    assert(apu->apu_output.sample_callback);
    apu->apu_output.sample_callback(apu, &filtered_output);
}

static uint16_t new_sweep_sample_legnth(struct apu_s *apu)
{
    uint16_t delta = apu->shadow_sweep_sample_legnth >> (apu->apu_registers[GB_IO_NR10] & 7);
    if (apu->apu_registers[GB_IO_NR10] & 8) {
        return apu->shadow_sweep_sample_legnth - delta;
    }
    return apu->shadow_sweep_sample_legnth + delta;
}

static void update_square_sample(struct apu_s *apu, unsigned index)
{
    if (apu->square_channels[index].current_sample_index & 0x80) return;

    uint8_t duty = apu->apu_registers[index == GB_SQUARE_1? GB_IO_NR11 :GB_IO_NR21] >> 6;
    update_sample(apu, index,
                  duties[apu->square_channels[index].current_sample_index + duty * 8]?
                  apu->square_channels[index].current_volume : 0,
                  0);
}


/* the effects of NRX2 writes on current volume are not well documented and differ
   between models and variants. The exact behavior can only be verified on CGB as it
   requires the PCM12 register. The behavior implemented here was verified on *my*
   CGB, which might behave differently from other CGB revisions, as well as from the
   DMG, MGB or SGB/2 */
static void nrx2_glitch(uint8_t *volume, uint8_t value, uint8_t old_value)
{
#if 0
    if (value & 8) {
        (*volume)++;
    }

    if (((value ^ old_value) & 8)) {
        *volume = 0x10 - *volume;
    }

    if ((value & 7) && !(old_value & 7) && *volume && !(value & 8)) {
        (*volume)--;
    }

    if ((old_value & 7) && (value & 8)) {
        (*volume)--;
    }

    (*volume) &= 0xF;
#endif
}

static void tick_square_envelope(struct apu_s *apu, enum GB_CHANNELS index)
{
    uint8_t nrx2 = apu->apu_registers[index == GB_SQUARE_1? GB_IO_NR12 : GB_IO_NR22];

    if (apu->square_channels[index].volume_countdown || (nrx2 & 7)) {
        if (!apu->square_channels[index].volume_countdown || !--apu->square_channels[index].volume_countdown) {
            if ((nrx2 & 8) && apu->square_channels[index].current_volume < 0xF) {
                apu->square_channels[index].current_volume++;
            }

            else if (!(nrx2 & 8) && apu->square_channels[index].current_volume > 0) {
                apu->square_channels[index].current_volume--;
            }

            apu->square_channels[index].volume_countdown = nrx2 & 7;

            if (apu->is_active[index]) {
                update_square_sample(apu, index);
            }
        }
    }
}

static void tick_noise_envelope(struct apu_s *apu)
{
    uint8_t nr42 = apu->apu_registers[GB_IO_NR42];

    if (apu->noise_channel.volume_countdown || (nr42 & 7)) {
        if (!--apu->noise_channel.volume_countdown) {
            if ((nr42 & 8) && apu->noise_channel.current_volume < 0xF) {
                apu->noise_channel.current_volume++;
            }

            else if (!(nr42 & 8) && apu->noise_channel.current_volume > 0) {
                apu->noise_channel.current_volume--;
            }

            apu->noise_channel.volume_countdown = nr42 & 7;

            if (apu->is_active[GB_NOISE]) {
                update_sample(apu, GB_NOISE,
                              (apu->noise_channel.lfsr & 1) ?
                              apu->noise_channel.current_volume : 0,
                              0);
            }
        }
    }
}

void GB_apu_div_event(struct apu_s *apu)
{
    if (!apu->global_enable) return;
    if (apu->skip_div_event) {
        apu->skip_div_event = false;
        return;
    }
    apu->div_divider++;

    if ((apu->div_divider & 1) == 0) {
        for (unsigned i = GB_SQUARE_2 + 1; i--;) {
            uint8_t nrx2 = apu->apu_registers[i == GB_SQUARE_1? GB_IO_NR12 : GB_IO_NR22];
            if (apu->is_active[i] && apu->square_channels[i].volume_countdown == 0 && (nrx2 & 7)) {
                tick_square_envelope(apu, i);
            }
        }

        if (apu->is_active[GB_NOISE] && apu->noise_channel.volume_countdown == 0 && (apu->apu_registers[GB_IO_NR42] & 7)) {
            tick_noise_envelope(apu);
        }
    }

    if ((apu->div_divider & 7) == 0) {
        for (unsigned i = GB_SQUARE_2 + 1; i--;) {
            tick_square_envelope(apu, i);
        }

        tick_noise_envelope(apu);
    }

    if ((apu->div_divider & 1) == 1) {
        for (unsigned i = GB_SQUARE_2 + 1; i--;) {
            if (apu->square_channels[i].length_enabled) {
                if (apu->square_channels[i].pulse_length) {
                    if (!--apu->square_channels[i].pulse_length) {
                        apu->is_active[i] = false;
                        update_sample(apu, i, 0, 0);
                    }
                }
            }
        }

        if (apu->wave_channel.length_enabled) {
            if (apu->wave_channel.pulse_length) {
                if (!--apu->wave_channel.pulse_length) {
                    apu->is_active[GB_WAVE] = false;
                    update_sample(apu, GB_WAVE, 0, 0);
                }
            }
        }

        if (apu->noise_channel.length_enabled) {
            if (apu->noise_channel.pulse_length) {
                if (!--apu->noise_channel.pulse_length) {
                    apu->is_active[GB_NOISE] = false;
                    update_sample(apu, GB_NOISE, 0, 0);
                }
            }
        }
    }

    if ((apu->div_divider & 3) == 3) {
        if (!apu->sweep_enabled) {
            return;
        }
        if (apu->square_sweep_countdown) {
            if (!--apu->square_sweep_countdown) {
                if ((apu->apu_registers[GB_IO_NR10] & 0x70) && (apu->apu_registers[GB_IO_NR10] & 0x07)) {
                    apu->square_channels[GB_SQUARE_1].sample_length =
                        apu->shadow_sweep_sample_legnth =
                        apu->new_sweep_sample_legnth;
                }

                if (apu->apu_registers[GB_IO_NR10] & 0x70) {
                    /* Recalculation and overflow check only occurs after a delay */
                    apu->square_sweep_calculate_countdown = 0x13 - apu->lf_div;
                }

                apu->square_sweep_countdown = ((apu->apu_registers[GB_IO_NR10] >> 4) & 7);
                if (!apu->square_sweep_countdown) apu->square_sweep_countdown = 8;
            }
        }
    }
}

void GB_apu_run(struct apu_s *apu)
{
	/* Convert 4MHZ to 2MHz. apu_cycles is always divisable by 4. */
	uint8_t apu_cycles = apu->apu_cycles >> 2;
	apu->apu_cycles = 0;
	apu->apu_output.sample_cycles += apu_cycles;

	/* TODO: Probably not required. */
	if (!apu_cycles)
		return;

	/* To align the square signal to 1MHz */
	apu->lf_div ^= apu_cycles & 1;
	apu->noise_channel.alignment += apu_cycles;

	if (apu->square_sweep_calculate_countdown) {
		if (apu->square_sweep_calculate_countdown > apu_cycles) {
			apu->square_sweep_calculate_countdown -= apu_cycles;
		}
		else {
			/* APU bug: sweep frequency is checked after adding the sweep delta twice */
			apu->new_sweep_sample_legnth = new_sweep_sample_legnth(apu);
			if (apu->new_sweep_sample_legnth > 0x7ff) {
				apu->is_active[GB_SQUARE_1] = false;
				update_sample(apu, GB_SQUARE_1, 0, apu->square_sweep_calculate_countdown - apu_cycles);
				apu->sweep_enabled = false;
			}
			apu->sweep_decreasing |= apu->apu_registers[GB_IO_NR10] & 8;
			apu->square_sweep_calculate_countdown = 0;
		}
	}

	UNROLL
	for (unsigned i = GB_SQUARE_1; i <= GB_SQUARE_2; i++) {
		if (apu->is_active[i]) {
			uint8_t apu_cycles_left = apu_cycles;
			while (unlikely(apu_cycles_left > apu->square_channels[i].sample_countdown)) {
				apu_cycles_left -= apu->square_channels[i].sample_countdown + 1;
				apu->square_channels[i].sample_countdown = (apu->square_channels[i].sample_length ^ 0x7FF) * 2 + 1;
				apu->square_channels[i].current_sample_index++;
				apu->square_channels[i].current_sample_index &= 0x7;

				update_square_sample(apu, i);
			}
			if (apu_cycles_left) {
				apu->square_channels[i].sample_countdown -= apu_cycles_left;
			}
		}
	}

	apu->wave_channel.wave_form_just_read = false;
	if (apu->is_active[GB_WAVE]) {
		uint8_t apu_cycles_left = apu_cycles;
		while (unlikely(apu_cycles_left > apu->wave_channel.sample_countdown)) {
			apu_cycles_left -= apu->wave_channel.sample_countdown + 1;
			apu->wave_channel.sample_countdown = apu->wave_channel.sample_length ^ 0x7FF;
			apu->wave_channel.current_sample_index++;
			apu->wave_channel.current_sample_index &= 0x1F;
			apu->wave_channel.current_sample =
				apu->wave_channel.wave_form[apu->wave_channel.current_sample_index];
			update_sample(apu, GB_WAVE,
					apu->wave_channel.current_sample >> apu->wave_channel.shift,
					apu_cycles - apu_cycles_left);
			apu->wave_channel.wave_form_just_read = true;
		}
		if (apu_cycles_left) {
			apu->wave_channel.sample_countdown -= apu_cycles_left;
			apu->wave_channel.wave_form_just_read = false;
		}
	}

	if (apu->is_active[GB_NOISE]) {
		uint8_t apu_cycles_left = apu_cycles;
		while (unlikely(apu_cycles_left > apu->noise_channel.sample_countdown)) {
			apu_cycles_left -= apu->noise_channel.sample_countdown + 1;
			apu->noise_channel.sample_countdown = apu->noise_channel.sample_length * 4 + 3;

			/* Step LFSR */
			unsigned high_bit_mask = apu->noise_channel.narrow ? 0x4040 : 0x4000;
			/* Todo: is this formula is different on a GBA? */
			bool new_high_bit = (apu->noise_channel.lfsr ^ (apu->noise_channel.lfsr >> 1) ^ 1) & 1;
			apu->noise_channel.lfsr >>= 1;

			if (new_high_bit) {
				apu->noise_channel.lfsr |= high_bit_mask;
			}
			else {
				/* This code is not redundent, it's relevant when switching LFSR widths */
				apu->noise_channel.lfsr &= ~high_bit_mask;
			}

			apu->current_lfsr_sample = apu->noise_channel.lfsr & 1;

			update_sample(apu, GB_NOISE,
					apu->current_lfsr_sample ?
					apu->noise_channel.current_volume : 0,
					0);
		}
		if (apu_cycles_left) {
			apu->noise_channel.sample_countdown -= apu_cycles_left;
		}
	}

	if (apu->apu_output.sample_rate) {
		apu->apu_output.cycles_since_render += apu_cycles;

		// TODO: Shouldn't this be a for loop?
		if (apu->apu_output.sample_cycles >= apu->apu_output.cycles_per_sample) {
			apu->apu_output.sample_cycles -= apu->apu_output.cycles_per_sample;
			render(apu);
		}
	}
}
void GB_apu_init(struct apu_s *apu)
{
	//memset(apu, 0, sizeof(*apu));
	/* Restore the wave form */
	for (unsigned reg = GB_IO_WAV_START; reg <= GB_IO_WAV_END; reg++)
	{
		apu->wave_channel.wave_form[(reg - GB_IO_WAV_START) * 2]     = apu->apu_registers[reg] >> 4;
		apu->wave_channel.wave_form[(reg - GB_IO_WAV_START) * 2 + 1] = apu->apu_registers[reg] & 0xF;
	}

	apu->lf_div = 1;
#if 0
	/* APU glitch: When turning the APU on while DIV's bit 4 (or 5 in double speed mode) is on,
	   the first DIV/APU event is skipped. */
	if (apu->div_counter & (apu->cgb_double_speed? 0x2000 : 0x1000)) {
		apu->skip_div_event = true;
	}
#endif
}

uint8_t GB_apu_read(struct apu_s *apu, uint8_t reg)
{
    if (reg == GB_IO_NR52) {
        uint8_t value = 0;
        for (int i = 0; i < GB_N_CHANNELS; i++) {
            value >>= 1;
            if (apu->is_active[i]) {
                value |= 0x8;
            }
        }
        if (apu->global_enable) {
            value |= 0x80;
        }
        value |= 0x70;
        return value;
    }

    static const char read_mask[GB_IO_WAV_END - GB_IO_NR10 + 1] = {
     /* NRX0  NRX1  NRX2  NRX3  NRX4 */
        0x80, 0x3F, 0x00, 0xFF, 0xBF, // NR1X
        0xFF, 0x3F, 0x00, 0xFF, 0xBF, // NR2X
        0x7F, 0xFF, 0x9F, 0xFF, 0xBF, // NR3X
        0xFF, 0xFF, 0x00, 0x00, 0xBF, // NR4X
        0x00, 0x00, 0x70, 0xFF, 0xFF, // NR5X

        0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, // Unused
        // Wave RAM
        0, /* ... */
    };

    if (reg >= GB_IO_WAV_START && reg <= GB_IO_WAV_END && apu->is_active[GB_WAVE]) {
#if 0
        if (!GB_is_cgb(apu) && !apu->wave_channel.wave_form_just_read) {
            return 0xFF;
        }
#endif
        reg = GB_IO_WAV_START + apu->wave_channel.current_sample_index / 2;
    }

    return apu->apu_registers[reg] | read_mask[reg - GB_IO_NR10];
}

void GB_apu_write(struct apu_s *apu, uint8_t reg, uint8_t value)
{
	if (!apu->global_enable &&
			reg != GB_IO_NR52 &&
			reg < GB_IO_WAV_START)
	{
		return;
	}

    if (reg >= GB_IO_WAV_START && reg <= GB_IO_WAV_END && apu->is_active[GB_WAVE]) {
#if 0
        if (!GB_is_cgb(apu) && !apu->wave_channel.wave_form_just_read) {
            return;
        }
#endif
        reg = GB_IO_WAV_START + apu->wave_channel.current_sample_index / 2;
    }

    /* Todo: this can and should be rewritten with a function table. */
    switch (reg) {
        /* Globals */
        case GB_IO_NR50:
        case GB_IO_NR51:
            apu->apu_registers[reg] = value;
            /* These registers affect the output of all 4 channels (but not the output of the PCM registers).*/
            /* We call update_samples with the current value so the APU output is updated with the new outputs */
            for (unsigned i = GB_N_CHANNELS; i--;) {
                update_sample(apu, i, apu->samples[i], 0);
            }
            break;
        case GB_IO_NR52: {

            uint8_t old_nrx1[] = {
                apu->apu_registers[GB_IO_NR11],
                apu->apu_registers[GB_IO_NR21],
                apu->apu_registers[GB_IO_NR31],
                apu->apu_registers[GB_IO_NR41]
            };
            if ((value & 0x80) && !apu->global_enable) {
                GB_apu_init(apu);
                apu->global_enable = true;
            }
            else if (!(value & 0x80) && apu->global_enable)  {
                for (unsigned i = GB_N_CHANNELS; i--;) {
                    update_sample(apu, i, 0, 0);
                }
                //memset(apu, 0, sizeof(*apu));
                memset(apu->apu_registers + GB_IO_NR10, 0, GB_IO_WAV_START - GB_IO_NR10);
                old_nrx1[0] &= 0x3F;
                old_nrx1[1] &= 0x3F;

                apu->global_enable = false;
            }

#if 0
            if (!GB_is_cgb(apu) && (value & 0x80)) {
                GB_apu_write(apu, GB_IO_NR11, old_nrx1[0]);
                GB_apu_write(apu, GB_IO_NR21, old_nrx1[1]);
                GB_apu_write(apu, GB_IO_NR31, old_nrx1[2]);
                GB_apu_write(apu, GB_IO_NR41, old_nrx1[3]);
            }
#endif
        }
        break;

        /* Square channels */
        case GB_IO_NR10:
            if (apu->sweep_decreasing && !(value & 8)) {
                apu->is_active[GB_SQUARE_1] = false;
                update_sample(apu, GB_SQUARE_1, 0, 0);
                apu->sweep_enabled = false;
                apu->square_sweep_calculate_countdown = 0;
            }
            if ((value & 0x70) == 0) {
                /* Todo: what happens if we set period to 0 while a calculate event is scheduled, and then
                         re-set it to non-zero? */
                apu->square_sweep_calculate_countdown = 0;
            }
            break;

        case GB_IO_NR11:
        case GB_IO_NR21: {
            unsigned index = reg == GB_IO_NR21? GB_SQUARE_2: GB_SQUARE_1;
            apu->square_channels[index].pulse_length = (0x40 - (value & 0x3f));
            if (!apu->global_enable) {
                value &= 0x3f;
            }
            break;
        }

        case GB_IO_NR12:
        case GB_IO_NR22: {
            unsigned index = reg == GB_IO_NR22? GB_SQUARE_2: GB_SQUARE_1;
            if (((value & 0x7) == 0) && ((apu->apu_registers[reg] & 0x7) != 0)) {
                /* Envelope disabled */
                apu->square_channels[index].volume_countdown = 0;
            }
            if ((value & 0xF8) == 0) {
                /* This disables the DAC */
                apu->apu_registers[reg] = value;
                apu->is_active[index] = false;
                update_sample(apu, index, 0, 0);
            }
            else if (apu->is_active[index]) {
                nrx2_glitch(&apu->square_channels[index].current_volume, value, apu->apu_registers[reg]);
                update_square_sample(apu, index);
            }

            break;
        }

        case GB_IO_NR13:
        case GB_IO_NR23: {
            unsigned index = reg == GB_IO_NR23? GB_SQUARE_2: GB_SQUARE_1;
            apu->square_channels[index].sample_length &= ~0xFF;
            apu->square_channels[index].sample_length |= value & 0xFF;
            break;
        }

        case GB_IO_NR14:
        case GB_IO_NR24: {
            unsigned index = reg == GB_IO_NR24? GB_SQUARE_2: GB_SQUARE_1;
            apu->square_channels[index].sample_length &= 0xFF;
            apu->square_channels[index].sample_length |= (value & 7) << 8;
            if (index == GB_SQUARE_1) {
                apu->shadow_sweep_sample_legnth =
                    apu->new_sweep_sample_legnth =
                    apu->square_channels[0].sample_length;
            }
            if (value & 0x80) {
                /* Current sample index remains unchanged when restarting channels 1 or 2. It is only reset by
                   turning the APU off. */
                if (!apu->is_active[index]) {
                    apu->square_channels[index].sample_countdown = (apu->square_channels[index].sample_length ^ 0x7FF) * 2 + 6 - apu->lf_div;
                }
                else {
                    /* Timing quirk: if already active, sound starts 2 (2MHz) ticks earlier.*/
                    apu->square_channels[index].sample_countdown = (apu->square_channels[index].sample_length ^ 0x7FF) * 2 + 4 - apu->lf_div;
                }
                apu->square_channels[index].current_volume = apu->apu_registers[index == GB_SQUARE_1 ? GB_IO_NR12 : GB_IO_NR22] >> 4;

                /* The volume changes caused by NRX4 sound start take effect instantly (i.e. the effect the previously
                   started sound). The playback itself is not instant which is why we don't update the sample for other
                   cases. */
                if (apu->is_active[index]) {
                    update_square_sample(apu, index);
                }

                apu->square_channels[index].volume_countdown = apu->apu_registers[index == GB_SQUARE_1 ? GB_IO_NR12 : GB_IO_NR22] & 7;

                if ((apu->apu_registers[index == GB_SQUARE_1 ? GB_IO_NR12 : GB_IO_NR22] & 0xF8) != 0 && !apu->is_active[index]) {
                    apu->is_active[index] = true;
                    update_sample(apu, index, 0, 0);
                    /* We use the highest bit in current_sample_index to mark this sample is not actually playing yet, */
                    apu->square_channels[index].current_sample_index |= 0x80;
                }
                if (apu->square_channels[index].pulse_length == 0) {
                    apu->square_channels[index].pulse_length = 0x40;
                    apu->square_channels[index].length_enabled = false;
                }

                if (index == GB_SQUARE_1) {
                    apu->sweep_decreasing = false;
                    if (apu->apu_registers[GB_IO_NR10] & 7) {
                        /* APU bug: if shift is nonzero, overflow check also occurs on trigger */
                        apu->square_sweep_calculate_countdown = 0x13 - apu->lf_div;
                    }
                    else {
                        apu->square_sweep_calculate_countdown = 0;
                    }
                    apu->sweep_enabled = apu->apu_registers[GB_IO_NR10] & 0x77;
                    apu->square_sweep_countdown = ((apu->apu_registers[GB_IO_NR10] >> 4) & 7);
                    if (!apu->square_sweep_countdown) apu->square_sweep_countdown = 8;
                }

            }

            /* APU glitch - if length is enabled while the DIV-divider's LSB is 1, tick the length once. */
            if ((value & 0x40) &&
                !apu->square_channels[index].length_enabled &&
                (apu->div_divider & 1) &&
                apu->square_channels[index].pulse_length) {
                apu->square_channels[index].pulse_length--;
                if (apu->square_channels[index].pulse_length == 0) {
                    if (value & 0x80) {
                        apu->square_channels[index].pulse_length = 0x3F;
                    }
                    else {
                        apu->is_active[index] = false;
                        update_sample(apu, index, 0, 0);
                    }
                }
            }
            apu->square_channels[index].length_enabled = value & 0x40;
            break;
        }

        /* Wave channel */
        case GB_IO_NR30:
            apu->wave_channel.enable = value & 0x80;
            if (!apu->wave_channel.enable) {
                apu->is_active[GB_WAVE] = false;
                update_sample(apu, GB_WAVE, 0, 0);
            }
            break;
        case GB_IO_NR31:
            apu->wave_channel.pulse_length = (0x100 - value);
            break;
        case GB_IO_NR32:
            apu->wave_channel.shift = (uint8_t[]){4, 0, 1, 2}[(value >> 5) & 3];
            if (apu->is_active[GB_WAVE]) {
                update_sample(apu, GB_WAVE, apu->wave_channel.current_sample >> apu->wave_channel.shift, 0);
            }
            break;
        case GB_IO_NR33:
            apu->wave_channel.sample_length &= ~0xFF;
            apu->wave_channel.sample_length |= value & 0xFF;
            break;
	case GB_IO_NR34:
	    apu->wave_channel.sample_length &= 0xFF;
	    apu->wave_channel.sample_length |= (value & 7) << 8;
	    if ((value & 0x80)) {
#if 0
		    /* DMG bug: wave RAM gets corrupted if the channel is retriggerred 1 cycle before the APU
		       reads from it. */
		    if (!GB_is_cgb(apu) &&
				    apu->is_active[GB_WAVE] &&
				    apu->wave_channel.sample_countdown == 0 &&
				    apu->wave_channel.enable) {
			    unsigned offset = ((apu->wave_channel.current_sample_index + 1) >> 1) & 0xF;

			    /* This glitch varies between models and even specific instances:
			       DMG-B:     Most of them behave as emulated. A few behave differently.
SGB:       As far as I know, all tested instances behave as emulated.
MGB, SGB2: Most instances behave non-deterministically, a few behave as emulated.

Additionally, I believe DMGs, including those we behave differently than emulated,
are all deterministic. */
			    if (offset < 4) {
				    apu->apu_registers[GB_IO_WAV_START] = apu->apu_registers[GB_IO_WAV_START + offset];
				    apu->wave_channel.wave_form[0] = apu->wave_channel.wave_form[offset / 2];
				    apu->wave_channel.wave_form[1] = apu->wave_channel.wave_form[offset / 2 + 1];
			    }
			    else {
				    memcpy(apu->apu_registers + GB_IO_WAV_START,
						    apu->apu_registers + GB_IO_WAV_START + (offset & ~3),
						    4);
				    memcpy(apu->wave_channel.wave_form,
						    apu->wave_channel.wave_form + (offset & ~3) * 2,
						    8);
			    }
		    }
#endif
		    if (!apu->is_active[GB_WAVE]) {
			    apu->is_active[GB_WAVE] = true;
			    update_sample(apu, GB_WAVE,
					    apu->wave_channel.current_sample >> apu->wave_channel.shift,
					    0);
		    }
		    apu->wave_channel.sample_countdown = (apu->wave_channel.sample_length ^ 0x7FF) + 3;
		    apu->wave_channel.current_sample_index = 0;
		    if (apu->wave_channel.pulse_length == 0) {
			    apu->wave_channel.pulse_length = 0x100;
			    apu->wave_channel.length_enabled = false;
		    }
		    /* Note that we don't change the sample just yet! This was verified on hardware. */
	    }

	    /* APU glitch - if length is enabled while the DIV-divider's LSB is 1, tick the length once. */
	    if ((value & 0x40) &&
			    !apu->wave_channel.length_enabled &&
			    (apu->div_divider & 1) &&
			    apu->wave_channel.pulse_length) {
		    apu->wave_channel.pulse_length--;
		    if (apu->wave_channel.pulse_length == 0) {
			    if (value & 0x80) {
				    apu->wave_channel.pulse_length = 0xFF;
			    }
			    else {
				    apu->is_active[GB_WAVE] = false;
				    update_sample(apu, GB_WAVE, 0, 0);
			    }
		    }
	    }
	    apu->wave_channel.length_enabled = value & 0x40;
	    if (apu->is_active[GB_WAVE] && !apu->wave_channel.enable) {
		    apu->is_active[GB_WAVE] = false;
		    update_sample(apu, GB_WAVE, 0, 0);
	    }

	    break;

        /* Noise Channel */

        case GB_IO_NR41: {
            apu->noise_channel.pulse_length = (0x40 - (value & 0x3f));
            break;
        }

        case GB_IO_NR42: {
            if (((value & 0x7) == 0) && ((apu->apu_registers[reg] & 0x7) != 0)) {
                /* Envelope disabled */
                apu->noise_channel.volume_countdown = 0;
            }
            if ((value & 0xF8) == 0) {
                /* This disables the DAC */
                apu->apu_registers[reg] = value;
                apu->is_active[GB_NOISE] = false;
                update_sample(apu, GB_NOISE, 0, 0);
            }
            else if (apu->is_active[GB_NOISE]){
                nrx2_glitch(&apu->noise_channel.current_volume, value, apu->apu_registers[reg]);
                update_sample(apu, GB_NOISE,
                              apu->current_lfsr_sample ?
                              apu->noise_channel.current_volume : 0,
                              0);
            }
            break;
        }

        case GB_IO_NR43: {
            apu->noise_channel.narrow = value & 8;
            unsigned divisor = (value & 0x07) << 1;
            if (!divisor) divisor = 1;
            apu->noise_channel.sample_length = (divisor << (value >> 4)) - 1;

            /* Todo: changing the frequency sometimes delays the next sample. This is probably
               due to how the frequency is actually calculated in the noise channel, which is probably
               not by calculating the effective sample length and counting simiarly to the other channels.
               This is not emulated correctly. */
            break;
        }

        case GB_IO_NR44: {
            if (value & 0x80) {
                apu->noise_channel.sample_countdown = (apu->noise_channel.sample_length) * 2 + 6 - apu->lf_div;

                /* I'm COMPLETELY unsure about this logic, but it passes all relevant tests.
                   See comment in NR43. */
                if ((apu->apu_registers[GB_IO_NR43] & 7) && (apu->noise_channel.alignment & 2) == 0) {
                    if ((apu->apu_registers[GB_IO_NR43] & 7) == 1) {
                        apu->noise_channel.sample_countdown += 2;
                    }
                    else {
                        apu->noise_channel.sample_countdown -= 2;
                    }
                }
                if (apu->is_active[GB_NOISE]) {
                    apu->noise_channel.sample_countdown += 2;
                }

                apu->noise_channel.current_volume = apu->apu_registers[GB_IO_NR42] >> 4;

                /* The volume changes caused by NRX4 sound start take effect instantly (i.e. the effect the previously
                 started sound). The playback itself is not instant which is why we don't update the sample for other
                 cases. */
                if (apu->is_active[GB_NOISE]) {
                    update_sample(apu, GB_NOISE,
                                  apu->current_lfsr_sample ?
                                  apu->noise_channel.current_volume : 0,
                                  0);
                }
                apu->noise_channel.lfsr = 0;
                apu->current_lfsr_sample = false;
                apu->noise_channel.volume_countdown = apu->apu_registers[GB_IO_NR42] & 7;

                if (!apu->is_active[GB_NOISE] && (apu->apu_registers[GB_IO_NR42] & 0xF8) != 0) {
                    apu->is_active[GB_NOISE] = true;
                    update_sample(apu, GB_NOISE, 0, 0);
                }

                if (apu->noise_channel.pulse_length == 0) {
                    apu->noise_channel.pulse_length = 0x40;
                    apu->noise_channel.length_enabled = false;
                }
            }

            /* APU glitch - if length is enabled while the DIV-divider's LSB is 1, tick the length once. */
            if ((value & 0x40) &&
                !apu->noise_channel.length_enabled &&
                (apu->div_divider & 1) &&
                apu->noise_channel.pulse_length) {
                apu->noise_channel.pulse_length--;
                if (apu->noise_channel.pulse_length == 0) {
                    if (value & 0x80) {
                        apu->noise_channel.pulse_length = 0x3F;
                    }
                    else {
                        apu->is_active[GB_NOISE] = false;
                        update_sample(apu, GB_NOISE, 0, 0);
                    }
                }
            }
            apu->noise_channel.length_enabled = value & 0x40;
            break;
        }

        default:
            if (reg >= GB_IO_WAV_START && reg <= GB_IO_WAV_END) {
                apu->wave_channel.wave_form[(reg - GB_IO_WAV_START) * 2]     = value >> 4;
                apu->wave_channel.wave_form[(reg - GB_IO_WAV_START) * 2 + 1] = value & 0xF;
            }
    }
    apu->apu_registers[reg] = value;
}

static void GB_apu_update_cycles_per_sample(struct apu_s *apu)
{
    if (apu->apu_output.sample_rate) {
        apu->apu_output.cycles_per_sample =
		CPU_FREQUENCY / (double)apu->apu_output.sample_rate;
    }
}

void GB_set_sample_rate(struct apu_s *apu, const uint_fast16_t sample_rate)
{
    apu->apu_output.sample_rate = sample_rate;
#if 1
    if (sample_rate) {
        apu->apu_output.highpass_rate = pow(0.999958, CPU_FREQUENCY / (double)sample_rate);
    }
#endif
    GB_apu_update_cycles_per_sample(apu);
}

void GB_apu_set_sample_callback(struct apu_s *apu, GB_sample_callback_t callback)
{
    apu->apu_output.sample_callback = callback;
}

#if 0
void GB_set_highpass_filter_mode(struct apu_s *apu, GB_highpass_mode_t mode)
{
    apu->apu_output.highpass_mode = mode;
}
#endif
