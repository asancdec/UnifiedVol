# Architecture

The repository is organized as a layered library: examples and tests use the
public headers under `uv/`; model implementations depend on shared core,
mathematical, and optimization components; low-level callback or implementation
details are kept out of the public API.

In the diagram below, arrows mean "depends on, includes, or calls".

```mermaid
flowchart LR
    subgraph Clients["Clients"]
        direction TB
        Examples[examples]
        Tests[tests]
    end

    subgraph Public["Public API"]
        direction TB
        API[uv/UnifiedVol.hpp]
    end

    subgraph Domain["Domain Layer"]
        direction TB
        Models[uv/Models<br/>SVI, Heston]
        IO[uv/IO]
        Core[uv/Core]
    end

    subgraph Numerics["Numerical Layer"]
        direction TB
        Math[uv/Math]
        Opt[uv/Optimization]
    end

    subgraph Foundation["Foundation"]
        direction TB
        Base[uv/Base]
    end

    subgraph External["External"]
        direction TB
        ExtTest[GoogleTest]
        ExtOpt[NLopt / Ceres]
        ExtMath[Boost.Math<br/>lets_be_rational]
    end

    Examples --> API
    Tests --> API
    Tests --> ExtTest

    API --> Base
    API --> Core
    API --> IO
    API --> Models
    API --> Math
    API --> Opt

    Core --> Base

    IO --> Base
    IO --> Core
    IO --> Math

    Models --> Base
    Models --> Core
    Models --> Math
    Models --> Opt

    Math --> Base
    Math --> Core
    Math --> ExtMath

    Opt --> Base
    Opt --> ExtOpt

    Examples ~~~ Tests
    Models ~~~ IO
    Math ~~~ Opt
    Base ~~~ ExtTest
```
