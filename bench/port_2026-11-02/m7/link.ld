/* Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche). Apache License 2.0.
   Memory of QEMU's mps2-an500 (Cortex-M7): code in SSRAM1 at 0, data and stack in SSRAM2/3 at 0x20000000. */
MEMORY { FLASH (rx) : ORIGIN = 0x00000000, LENGTH = 4M  RAM (rwx) : ORIGIN = 0x20000000, LENGTH = 4M }
_estack = ORIGIN(RAM) + LENGTH(RAM);
SECTIONS {
  .isr_vector : { KEEP(*(.isr_vector)) } > FLASH
  .text : { *(.text*) *(.rodata*) KEEP(*(.init)) KEEP(*(.fini)) . = ALIGN(4);
            __preinit_array_start = .; KEEP(*(.preinit_array)) __preinit_array_end = .;
            __init_array_start = .; KEEP(*(.init_array*)) __init_array_end = .;
            __fini_array_start = .; KEEP(*(.fini_array*)) __fini_array_end = .; } > FLASH
  .ARM.exidx : { __exidx_start = .; *(.ARM.exidx*) __exidx_end = .; } > FLASH
  _sidata = LOADADDR(.data);
  .data : { _sdata = .; *(.data*) . = ALIGN(4); _edata = .; } > RAM AT > FLASH
  .bss (NOLOAD) : { _sbss = .; *(.bss*) *(COMMON) . = ALIGN(4); _ebss = .; } > RAM
  end = .; _end = .; __end__ = .;
}
