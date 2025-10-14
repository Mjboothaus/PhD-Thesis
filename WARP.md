# PhD-Thesis Project - Charged Fluids Near Interfaces

## Project Overview

This project is a reimplementation of PhD thesis research on "Charged fluids near interfaces: Integral equation theory" originally submitted in 1998. The current implementation uses Python and Streamlit to create an interactive web application for exploring fluid behaviour near interfaces using integral equation theory.

### Key Features

- Interactive Streamlit web application for fluid modelling
- Support for multiple fluid types (Lennard-Jones, Potassium Chloride, Liquid Water, 2-2 Electrolyte)
- Numerical solver using Newton-Krylov algorithms
- Visualisation of convergence and solution curves
- Bulk fluid correlation function calculations using PyOZ
- Memory usage monitoring and performance tracking

## Technology Stack

- **Language**: Python 3.13
- **Frontend**: Streamlit 1.37+
- **Scientific Computing**: NumPy, SciPy, Matplotlib, Plotly
- **Data Management**: Pandas, DuckDB
- **Build Tool**: UV (Python package manager)
- **Task Runner**: Just (justfile)
- **Containerisation**: Docker
- **Deployment**: Google Cloud Run (production at https://phd.databooth.com.au)

## Project Structure

```
PhD-Thesis/
├── src/                          # Main application source code
│   ├── Main.py                   # Streamlit main application entry point
│   ├── pages/                    # Streamlit application pages
│   │   ├── 0_💡_Theory.py        # Theoretical background page
│   │   ├── 1_🫙_Lennard_Jones.py # Lennard-Jones fluid calculations
│   │   ├── 2_🫙_Potassium_Chloride.py # KCl fluid calculations
│   │   ├── 3_💧_Liquid_Water.py  # Water fluid calculations
│   │   ├── 4_📋_2-2_Electrolyte.py # 2-2 electrolyte calculations
│   │   └── 9_📂_GitHub_(Source_code).py # Source code reference
│   ├── modelling.py              # Core mathematical modelling functions
│   ├── parameters.py             # Fluid parameter definitions and management
│   ├── numerics.py               # Numerical computation utilities
│   ├── plotting.py               # Plotting and visualisation functions
│   ├── sidebar.py                # Streamlit sidebar components
│   ├── bulk/                     # Bulk fluid calculations (PyOZ integration)
│   └── oo_refactoring/           # Object-oriented refactoring (work in progress)
├── data/                         # Data files and fluid parameters
│   ├── fluid_parameters.toml     # Fluid configuration parameters
│   ├── *.data                    # Correlation function data files
│   └── digitised CSV files       # Experimental data
├── docs/                         # Documentation and theory files
├── notebooks/                    # Jupyter notebooks for development and analysis
├── tests/                        # Test suite (pytest)
├── output/                       # Generated output files
├── justfile                      # Task automation (Just build tool)
├── pyproject.toml               # Python project configuration (UV)
├── Dockerfile                   # Container configuration
└── README.md                    # Project documentation
```

## Development Environment Setup

### Prerequisites

- Python 3.13
- UV package manager
- Just task runner
- Git

### Initial Setup

```bash
# Clone the repository
git clone [repository-url]
cd PhD-Thesis

# Initialise development environment
just init

# or manually sync with UV
uv sync
```

## Available Tasks (Just Commands)

### Development Environment

- `just init` - Initialise development environment with UV
- `just sync` - Synchronise environment with pyproject.toml
- `just update` - Update UV environment and freeze requirements
- `just install-dev` - Install development requirements

### Application

- `just app` - Run the main Streamlit application
- `just app-opt` - Run optimised version of the app
- `just update-st-config` - Update Streamlit configuration

### Bulk Fluid Calculations

- `just bulk-fluid-pyoz [input_file]` - Run PyOZ bulk fluid calculations

### Testing and Quality

- `just test` - Run all tests
- `just test-cov` - Run tests with coverage report
- `just type-check` - Run type checking with mypy
- `just lint` - Run code linting (black, flake8)
- `just format` - Format code (black, isort)
- `just quality` - Run all quality checks (format, type-check, lint, test)

### Docker and Deployment

- `just docker-build [project_name]` - Build Docker image
- `just docker-run [project_name] [port]` - Run Docker container
- `just container` - Build and run Docker container
- `just deploy-render` - Deploy to Render.com

### Utilities

- `just reqs` - Generate requirements files
- `just clean` - Clean up cache and temporary files

## Key Components

### Fluid Types Supported

1. **Lennard-Jones (LJ1, LJ2)** - Simple atomic fluids
2. **Potassium Chloride (KCl)** - Ionic fluid system
3. **Liquid Water (H2O)** - Complex molecular fluid
4. **2-2 Electrolyte** - Divalent ionic systems

### Core Mathematical Components

- **Integral Equation Theory** - Ornstein-Zernike equations with closure relations
- **Newton-Krylov Solver** - Nonlinear equation solving
- **Correlation Functions** - Radial distribution and direct correlation functions
- **Wall Interactions** - Fluid-surface interface calculations

### Data Flow

1. **Parameter Definition** - Fluid properties loaded from TOML configuration
2. **Bulk Calculations** - Correlation functions from PyOZ or stored data
3. **Solver Setup** - Initial conditions and numerical parameters
4. **Iterative Solution** - Newton-Krylov algorithm convergence
5. **Visualisation** - Interactive plots of results and convergence

## Development Guidelines

### Code Style

- Follow PEP 8 style guidelines
- Use type hints throughout the codebase
- Maintain Australian English in documentation and comments
- Use US spelling for variables, classes, and function names
- Format code using Black and isort

### Git Workflow

- Use standard git branching workflow for features and fixes
- Create pull requests for all changes
- Write complete yet succinct commit messages
- Follow conventional commit format where appropriate

### Testing

- Write tests for all new functionality
- Maintain test coverage above 80%
- Use pytest for testing framework
- Include integration tests for key workflows

## Known Issues and TODOs

### Current Limitations

1. **Charged Fluid Correlations** - Short-range correlations don't approach zero as expected
2. **OO Refactoring** - Object-oriented refactoring is only partially complete
3. **Bulk Integration** - Automatic bulk fluid code integration needs improvement
4. **Data Structures** - Discrete data structures could be optimised

### Development Priorities

1. Fix charged fluid correlation function behaviour
2. Complete object-oriented refactoring
3. Integrate bulk fluid calculations automatically
4. Optimise data structures for better performance
5. Implement automated solver output capture using Streamlit widgets

## Performance Considerations

- **Memory Monitoring** - Built-in memory usage tracking with psutil
- **Solver Optimisation** - Newton-Krylov algorithms for efficient convergence
- **Caching** - Correlation function data caching to avoid recalculation
- **Streamlit Optimisation** - Efficient plot rendering and state management

## Deployment

### Local Development

```bash
just app
# Application available at http://localhost:8501
```

### Docker Deployment

```bash
just container
# Application available at http://localhost:8080
```

### Production

- Deployed on Google Cloud Run
- Production URL: https://phd.databooth.com.au
- Uses containerised deployment with optimised configuration

## Dependencies

### Core Dependencies

- streamlit>=1.37.1
- scipy>=1.14.0
- numpy>=2.0.1
- pandas>=2.2.2
- matplotlib>=3.9.2
- plotly>=5.23.0
- psutil>=7.1.0 (memory monitoring)
- loguru>=0.7.3 (logging)

### Development Dependencies

- pytest>=8.0.0
- black>=24.1.1
- isort>=5.13.2
- flake8>=7.0.0
- mypy>=1.8.0

## License

MIT License

## References

- Original PhD thesis: "Charged fluids near interfaces: Integral equation theory" (1998)
- PyOZ project: http://pyoz.vrbka.net (bulk fluid calculations)
- Fork available: https://github.com/ctk3b/pyoz

## Maintenance Notes

This project is actively maintained and represents a modern reimplementation of classical integral equation theory for fluid systems. The codebase is designed to be educational and research-oriented, with emphasis on clear mathematical implementation and interactive exploration of theoretical concepts.

Regular updates focus on:
- Numerical accuracy improvements
- Performance optimisation
- User interface enhancements
- Documentation and educational content
- Code quality and maintainability

---

**Created with UV and Just - Cross-platform compatible development environment**