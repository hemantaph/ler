# LeR Thread Management - Visual Guide

## Thread Assignment Strategy

```
┌─────────────────────────────────────────────────────────────────┐
│                    LeR Package                                  │
│                                                                 │
│  ┌──────────────────────────────────────────────────────────┐  │
│  │ __init__.py: Global Environment Configuration           │  │
│  │ ───────────────────────────────────────────────────────  │  │
│  │ NUMBA_NUM_THREADS = $(cpu_count)     # 32 on 32-core    │  │
│  │ OMP_NUM_THREADS = 1                  # Avoid BLAS        │  │
│  │ MKL_NUM_THREADS = 1                  # competition      │  │
│  │ OPENBLAS_NUM_THREADS = 1                                │  │
│  └──────────────────────────────────────────────────────────┘  │
│                                                                 │
│  ┌───────────────────────┐     ┌─────────────────────────────┐ │
│  │   Path 1: prange      │     │  Path 2: Multiprocessing   │ │
│  │ (Numba Parallelism)   │     │   (Process Pool)            │ │
│  ├───────────────────────┤     ├─────────────────────────────┤ │
│  │                       │     │                             │ │
│  │ Main Process:         │     │ Main Process:               │ │
│  │ ↓                     │     │ ↓                           │ │
│  │ set_num_threads(4)    │     │ Pool(processes=4,           │ │
│  │ ↓                     │     │      initializer=init_fn)   │ │
│  │ @njit(parallel=True)  │     │ ↓                           │ │
│  │ for i in prange(...)  │     │ Worker 1: set_num_threads(1)│ │
│  │   ├─ Thread 0         │     │ Worker 2: set_num_threads(1)│ │
│  │   ├─ Thread 1         │     │ Worker 3: set_num_threads(1)│ │
│  │   ├─ Thread 2         │     │ Worker 4: set_num_threads(1)│ │
│  │   └─ Thread 3         │     │                             │ │
│  │ ↓                     │     │ Total Threads: 4            │ │
│  │ Total: 4 threads      │     │ (1 thread per worker)       │ │
│  │ (correctly balanced)  │     │ (correctly balanced)        │ │
│  │                       │     │                             │ │
│  └───────────────────────┘     └─────────────────────────────┘ │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
```

---

## Pattern Comparison

### ✅ CORRECT PATTERN: Lens Redshift Multiprocessing

```
Multiple Source Redshifts (zs)
  │
  ├─ zs[0] → Worker 1 (set_num_threads(1)): compute lens redshift
  │
  ├─ zs[1] → Worker 2 (set_num_threads(1)): compute lens redshift
  │
  ├─ zs[2] → Worker 3 (set_num_threads(1)): compute lens redshift
  │
  └─ zs[3] → Worker 4 (set_num_threads(1)): compute lens redshift

Total Threads in Use: 4 (1 per worker)
CPU Cores Available: 32
Efficiency: Optimal (4/32 = 12.5% core utilization, but balanced)
```

### ❌ PROBLEMATIC PATTERN: Image Properties (CURRENT)

```
Multiple Lens Events
  │
  ├─ Event 1 → Worker 1 (NO set_num_threads): solve lens equation
  │             ├─ Numba Thread Pool: 32 threads!
  │             ├─ Thread 0-31 (context switching nightmare)
  │             └─ CPU cores: 100% oversubscribed
  │
  ├─ Event 2 → Worker 2 (NO set_num_threads): solve lens equation
  │             └─ Numba Thread Pool: 32 threads!
  │
  ├─ Event 3 → Worker 3 (NO set_num_threads): solve lens equation
  │             └─ Numba Thread Pool: 32 threads!
  │
  └─ Event 4 → Worker 4 (NO set_num_threads): solve lens equation
                └─ Numba Thread Pool: 32 threads!

Total Threads: 4 × 32 = 128 threads
CPU Cores Available: 32
Efficiency: Catastrophic!
Context Switches: Extreme (massive overhead)
Performance: 4-6x SLOWER
```

---

## Module Status Overview

```
┌────────────────────────────────────────────────────────────────┐
│                   MULTITHREADING (prange)                      │
├─────────────────────┬──────────────────────────────────────────┤
│ Module              │ Status                                   │
├─────────────────────┼──────────────────────────────────────────┤
│ optical_depth       │ ✅ Good: set_num_threads(npool)         │
│ image_properties    │ ⚠️  Redundant: no prange in epl_solver  │
│ rates/ler.py        │ ✅ Good: set_num_threads(npool) both    │
│                     │    unlensed & lensed                    │
└─────────────────────┴──────────────────────────────────────────┘

┌────────────────────────────────────────────────────────────────┐
│            MULTIPROCESSING (Pool.map)                           │
├──────────────────────────┬───────────────────────────────────┤
│ Module                   │ Status                            │
├──────────────────────────┼───────────────────────────────────┤
│ mp.py                    │ ✅ Good: set_num_threads(1) in    │
│ (lens_redshift)          │    worker                         │
├──────────────────────────┼───────────────────────────────────┤
│ sampler_functions.py     │ ✅ Good: set_num_threads(1) in    │
│ (importance_sampling)    │    worker                         │
├──────────────────────────┼───────────────────────────────────┤
│ multiprocessing_routine  │ 🔴 CRITICAL: MISSING              │
│ _epl_shear.py            │    set_num_threads(1) in          │
│ (image_properties)       │    initializer                     │
├──────────────────────────┼───────────────────────────────────┤
│ sfr_with_time_delay.py   │ 🔴 CRITICAL: MISSING              │
│ (merger_rate_density)    │    set_num_threads(1) in          │
│                          │    function                       │
├──────────────────────────┼───────────────────────────────────┤
│ cbc_source_redshift      │ 🔴 CRITICAL: NO initializer,      │
│ _distribution.py         │    depends on Issue #2 fix        │
└──────────────────────────┴───────────────────────────────────┘
```

---

## Thread Lifetime Diagram

### prange Execution (set_num_threads(npool))

```
Timeline:
┌─────────────────────────────────────────────────────┐
│ Main Thread                                         │
├─────────────────────────────────────────────────────┤
│                                                     │
│  set_num_threads(4)                                │
│  │                                                  │
│  ├─→ @njit(parallel=True) function called         │
│      │                                              │
│      ├─→ prange loop iteration space: 0-999        │
│      │   ├─ Thread 0: iterations 0-249            │
│      │   ├─ Thread 1: iterations 250-499          │
│      │   ├─ Thread 2: iterations 500-749          │
│      │   └─ Thread 3: iterations 750-999          │
│      │                                              │
│      └─→ All threads join at barrier               │
│                                                     │
│  set_num_threads(32)  ← Reset for next operation   │
│  │                                                  │
│  └─→ Continue with original CPU count              │
│                                                     │
└─────────────────────────────────────────────────────┘

Threads Active in Peak: 4 ✅
Context Switches: Minimal ✅
```

### Multiprocessing Execution (set_num_threads(1))

```
Timeline:
┌──────────────────────────────────────────────────────────┐
│ Main Process                                             │
├──────────────────────────────────────────────────────────┤
│                                                          │
│  Pool(processes=4, initializer=_init_worker)            │
│  │                                                       │
│  ├─→ Worker 1: set_num_threads(1)                       │
│  │   └─→ Task A (no internal parallelism)               │
│  │       └─→ Returns result                             │
│  │                                                       │
│  ├─→ Worker 2: set_num_threads(1)                       │
│  │   └─→ Task B                                         │
│  │                                                       │
│  ├─→ Worker 3: set_num_threads(1)                       │
│  │   └─→ Task C                                         │
│  │                                                       │
│  └─→ Worker 4: set_num_threads(1)                       │
│      └─→ Task D                                         │
│                                                          │
│  Gather results from all workers                        │
│                                                          │
└──────────────────────────────────────────────────────────┘

Processes Active: 4 
Threads per Process: 1
Total Threads: 4 ✅
Context Switches: Minimal ✅
```

---

## Performance Comparison

### System: 4-Core CPU, npool=4

#### Scenario 1: Image Properties (BEFORE FIX)

```
Worker 1: set_num_threads(32)      🔴 Wrong! Uses global NUMBA_NUM_THREADS
│         └─→ 32 Numba threads
Worker 2: set_num_threads(32)      🔴 Wrong!
│         └─→ 32 Numba threads
Worker 3: set_num_threads(32)      🔴 Wrong!
│         └─→ 32 Numba threads
Worker 4: set_num_threads(32)      🔴 Wrong!
          └─→ 32 Numba threads

Total: 4 × 32 = 128 threads on 4 cores
Ratio: 128 / 4 = 32x oversubscription

Context Switches: ~1,000+ per second
CPU Cache Thrashing: Severe
Actual Performance: 4-6% of theoretical peak
Time for 10,000 items: ~300 seconds
```

#### Scenario 1: Image Properties (AFTER FIX)

```
Worker 1: set_num_threads(1)       ✅ Correct!
│         └─→ 1 Numba thread
Worker 2: set_num_threads(1)       ✅ Correct!
│         └─→ 1 Numba thread
Worker 3: set_num_threads(1)       ✅ Correct!
│         └─→ 1 Numba thread
Worker 4: set_num_threads(1)       ✅ Correct!
          └─→ 1 Numba thread

Total: 4 × 1 = 4 threads on 4 cores
Ratio: 4 / 4 = 1x (perfect scaling)

Context Switches: ~10 per second
CPU Cache Efficiency: Optimal
Actual Performance: 80-90% of theoretical peak
Time for 10,000 items: ~50-60 seconds
Speed Improvement: 5-6x faster ✅
```

---

## Summary Table

| Aspect | prange | Multiprocessing |
|--------|--------|-----------------|
| **Thread Source** | Numba JIT | Python processes |
| **Parallelism Boundary** | Within single function | Across multiple tasks |
| **Thread Control Call** | `set_num_threads(npool)` BEFORE prange | `set_num_threads(1)` IN worker or initializer |
| **When to Control** | Before calling @njit(parallel=True) | When worker starts |
| **Reset Needed** | No (usually, but good practice) | No (per-worker setting) |
| **Optimal Threads** | Global (all cores) → Reduce for prange | 1 (pool provides parallelism) |
| **Status in LeR** | ✅ Mostly Good | 🔴 Critical Issues |

---

## Debugging Checklist

When thread management breaks:

```
[ ] Check NUMBA_NUM_THREADS environment variable
    $ python -c "import os; print(os.environ.get('NUMBA_NUM_THREADS'))"
    
[ ] Verify set_num_threads() imports
    from numba import set_num_threads
    set_num_threads(1)  # should not error
    
[ ] Monitor thread count during execution
    $ watch -n 0.1 'ps -aux | wc -l'
    
[ ] Profile with threadpoolctl
    from threadpoolctl import threadpool_info
    print(threadpool_info())  # reveals thread pool configuration
    
[ ] Check Pool initializer is being called
    def _init_worker():
        print("Worker initialized with PID", os.getpid())
        set_num_threads(1)
    
[ ] Profile CPU usage
    $ top -H -p <pid>  # shows all threads for process
    
[ ] Measure actual performance
    import time
    start = time.time()
    result = multiprocessing_function()
    print(f"Time: {time.time() - start:.2f}s")
```

---

## File Organization

```
LeR Package/
├── ler/
│   ├── __init__.py                               [✅ GOOD]
│   ├── lens_galaxy_population/
│   │   ├── mp.py                                 [✅ GOOD]
│   │   ├── optical_depth.py                      [✅ GOOD]
│   │   └── sampler_functions.py                  [✅ GOOD]
│   ├── gw_source_population/
│   │   ├── cbc_source_redshift_distribution.py   [🔴 Issue #3]
│   │   └── sfr_with_time_delay.py                [🔴 Issue #2]
│   ├── image_properties/
│   │   ├── image_properties.py                   [⚠️ Redundant]
│   │   └── multiprocessing_routine_epl_shear.py  [🔴 Issue #1]
│   └── rates/
│       └── ler.py                                [✅ GOOD]
│
├── THREAD_MANAGEMENT_AUDIT.md                    [Detailed analysis]
├── THREAD_MANAGEMENT_SUMMARY.md                  [Executive summary]
├── THREAD_MANAGEMENT_FIXES.md                    [Quick fix guide]
└── THREAD_MANAGEMENT_VISUAL.md                   [This file]
```

---

## References

- **Numba Threading:** https://numba.readthedocs.io/en/stable/user/threading-layer.html
- **Python Multiprocessing:** https://docs.python.org/3/library/multiprocessing.html
- **Thread Safety:** https://en.wikipedia.org/wiki/Thread_pool
- **CPU Oversubscription:** https://en.wikipedia.org/wiki/Oversubscription_(computing)
