#!/bin/bash

# Orv-Mol Installation Script
# ===========================

set -e  # Exit on any error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Print functions
print_header() {
    echo -e "${BLUE}======================================${NC}"
    echo -e "${BLUE}  Orv-Mol Installation Script${NC}"
    echo -e "${BLUE}======================================${NC}"
    echo
}

print_success() {
    echo -e "${GREEN}✓ $1${NC}"
}

print_error() {
    echo -e "${RED}✗ $1${NC}"
}

print_warning() {
    echo -e "${YELLOW}⚠ $1${NC}"
}

print_info() {
    echo -e "${BLUE}ℹ $1${NC}"
}

# Check if running on supported OS
check_os() {
    print_info "Checking operating system..."
    if [[ "$OSTYPE" == "linux-gnu"* ]]; then
        OS="linux"
        print_success "Linux detected"
    elif [[ "$OSTYPE" == "darwin"* ]]; then
        OS="macos"
        print_success "macOS detected"
    else
        print_error "Unsupported operating system: $OSTYPE"
        exit 1
    fi
}

# Check Python version
check_python() {
    print_info "Checking Python version..."
    if command -v python3 &> /dev/null; then
        PYTHON_VERSION=$(python3 -c 'import sys; print(".".join(map(str, sys.version_info[:2])))')
        PYTHON_MAJOR=$(echo $PYTHON_VERSION | cut -d. -f1)
        PYTHON_MINOR=$(echo $PYTHON_VERSION | cut -d. -f2)
        
        if [[ $PYTHON_MAJOR -eq 3 && $PYTHON_MINOR -ge 8 ]]; then
            print_success "Python $PYTHON_VERSION detected"
            PYTHON_CMD="python3"
        else
            print_error "Python 3.8+ required, found $PYTHON_VERSION"
            exit 1
        fi
    else
        print_error "Python 3 not found"
        exit 1
    fi
}

# Install system dependencies
install_system_deps() {
    print_info "Installing system dependencies..."
    
    if [[ "$OS" == "linux" ]]; then
        # Check if we can install packages
        if command -v apt-get &> /dev/null; then
            if ! dpkg -l | grep -q autodock-vina; then
                print_info "Installing AutoDock Vina..."
                sudo apt-get update
                sudo apt-get install -y autodock-vina
                print_success "AutoDock Vina installed"
            else
                print_success "AutoDock Vina already installed"
            fi
        elif command -v yum &> /dev/null; then
            print_warning "YUM detected. Please install autodock-vina manually"
        else
            print_warning "Unknown package manager. Please install autodock-vina manually"
        fi
    elif [[ "$OS" == "macos" ]]; then
        if command -v brew &> /dev/null; then
            if ! brew list autodock-vina &> /dev/null; then
                print_info "Installing AutoDock Vina via Homebrew..."
                brew install autodock-vina
                print_success "AutoDock Vina installed"
            else
                print_success "AutoDock Vina already installed"
            fi
        else
            print_warning "Homebrew not found. Please install autodock-vina manually"
        fi
    fi
}

# Create virtual environment
create_venv() {
    print_info "Setting up Python virtual environment..."
    
    if [[ ! -d "venv" ]]; then
        $PYTHON_CMD -m venv venv
        print_success "Virtual environment created"
    else
        print_info "Virtual environment already exists"
    fi
    
    # Activate virtual environment
    source venv/bin/activate
    print_success "Virtual environment activated"
}

# Install Python dependencies
install_python_deps() {
    print_info "Installing Python dependencies..."
    
    # Upgrade pip
    pip install --upgrade pip
    
    # Install requirements
    if [[ -f "requirements.txt" ]]; then
        pip install -r requirements.txt
        print_success "Python dependencies installed"
    else
        print_error "requirements.txt not found"
        exit 1
    fi
}

# Verify installation
verify_installation() {
    print_info "Verifying installation..."
    
    # Check Vina
    if command -v vina &> /dev/null; then
        print_success "AutoDock Vina: $(vina --version 2>&1 | head -1)"
    else
        print_error "AutoDock Vina not found in PATH"
    fi
    
    # Check Python modules
    python3 -c "
import sys
modules = ['rdkit', 'Bio', 'numpy', 'pandas', 'matplotlib', 'yaml', 'loguru']
for module in modules:
    try:
        __import__(module)
        print(f'✓ {module}')
    except ImportError:
        print(f'✗ {module}')
        sys.exit(1)
"
    
    print_success "All modules verified"
}

# Run tests
run_tests() {
    print_info "Running installation tests..."
    
    if [[ -f "tests/test_installation.py" ]]; then
        python3 tests/test_installation.py
        print_success "Installation tests passed"
    else
        print_warning "Installation tests not found"
    fi
}

# Main installation function
main() {
    print_header
    
    check_os
    check_python
    install_system_deps
    create_venv
    install_python_deps
    verify_installation
    run_tests
    
    echo
    print_success "Installation completed successfully!"
    echo
    print_info "To get started:"
    echo "  1. Activate the virtual environment: source venv/bin/activate"
    echo "  2. Run a test: python dock.py --help"
    echo "  3. Try an example: python dock.py data/input/protein.pdb data/input/ligand_small.pdb"
    echo
    print_info "For more information, see README.md"
}

# Check if script is run directly
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    main "$@"
fi 