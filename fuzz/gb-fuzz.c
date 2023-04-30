#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#define ENABLE_LCD 1
#define ENABLE_SOUND 0

#include "../peanut_gb.h"

// Ram was max 64k
#define RAM_SIZE 65536

uint8_t ram[RAM_SIZE];

struct priv_t
{
	/* Pointer to allocated memory holding GB file. */
	const uint8_t *rom;
	size_t rom_size;
};

/**
 * Returns a byte from the ROM file at the given address.
 */
uint8_t gb_rom_read(struct gb_s *gb, const uint_fast32_t addr)
{
    const struct priv_t * const p = gb->direct.priv;
    return p->rom[addr % p->rom_size];
}

/**
 * Returns a byte from the cartridge RAM at the given address.
 */
uint8_t gb_cart_ram_read(struct gb_s *gb, const uint_fast32_t addr)
{
	return ram[addr % RAM_SIZE];
}

/**
 * Writes a given byte to the cartridge RAM at the given address.
 */
void gb_cart_ram_write(struct gb_s *gb, const uint_fast32_t addr,
	const uint8_t val)
{
	ram[addr % RAM_SIZE] = val;
}

/**
 * Handles an error reported by the emulator. The emulator context may be used
 * to better understand why the error given in gb_err was reported.
 */
void gb_error(struct gb_s *gb, const enum gb_error_e gb_err, const uint16_t val)
{
    if (gb_err == GB_HALT_FOREVER) {
        // Force this to stop
        exit(EXIT_FAILURE);
    }
}

int LLVMFuzzerTestOneInput(const uint8_t *Data, size_t Size) {
    if (Size > 0) {
        struct gb_s gb;
        struct priv_t priv;
        priv.rom = Data;
        priv.rom_size = Size;
        enum gb_init_error_e ret = gb_init(&gb, &gb_rom_read, &gb_cart_ram_read, &gb_cart_ram_write,
                &gb_error, &priv);

        if(ret != GB_INIT_NO_ERROR)
        {
            return 0;
        }
        
        uint8_t counter = 0;
        gb.gb_frame = 0;
        while (counter < 24 && gb.gb_frame == 0) {
            __gb_step_cpu(&gb);
            counter++;
        }
    }

    return 0;
}