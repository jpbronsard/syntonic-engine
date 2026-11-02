/* Copyright 2026 Jean-Pierre Bronsard (SyntonicAI Recherche)
   Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE. */
/* minimal Cortex-M7 startup for QEMU mps2-an500: vector table, FPU on, .data/.bss, C++ constructors, main, exit */
#include <stdint.h>
extern uint32_t _estack, _sidata, _sdata, _edata, _sbss, _ebss;
extern int main(void);
extern void __libc_init_array(void);
extern void initialise_monitor_handles(void);
extern void exit(int);
void Reset_Handler(void) {
  *(volatile uint32_t*)0xE000ED88 |= (0xFu << 20);          /* CPACR: CP10, CP11 full access */
  __asm volatile("dsb\n isb");
  uint32_t *s = &_sidata, *d = &_sdata; while (d < &_edata) *d++ = *s++;
  for (d = &_sbss; d < &_ebss;) *d++ = 0;
  initialise_monitor_handles();
  __libc_init_array();
  exit(main());
}
void Default_Handler(void) { for (;;) {} }
__attribute__((section(".isr_vector"), used)) void (*const vectors[16])(void) = {
  (void (*)(void))&_estack, Reset_Handler, Default_Handler, Default_Handler, Default_Handler, Default_Handler,
  Default_Handler, 0, 0, 0, 0, Default_Handler, Default_Handler, 0, Default_Handler, Default_Handler};
void _init(void) {} void _fini(void) {}
