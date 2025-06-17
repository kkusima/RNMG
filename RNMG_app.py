# Aurthor: Kenneth Kusima
# Version 1.0
# 06/2025
# =============================================================================
# REACTION NETWORK MATRIX GENERATOR
# =============================================================================
# This app helps create stoichiometric and atomic matrices from chemical 
# reaction mechanisms. It's designed for surface chemistry but works for 
# any type of reactions!

# Import all the libraries we need
import streamlit as st
import pandas as pd
import numpy as np
import re
from io import StringIO
from collections import defaultdict, Counter
import math

# =============================================================================
# PAGE SETUP AND CONFIGURATION
# =============================================================================

# Set up the page layout and basic configuration
st.set_page_config(
    page_title="Reaction Network Matrix Generator",
    page_icon="🔗",
    layout="wide"  # Uses full width of browser
)

# Main title and description
st.title("Reaction Network Matrix Generator")
st.markdown("Create stoichiometric and atomic matrices from chemical reaction mechanisms")

# =============================================================================
# SESSION STATE INITIALIZATION
# =============================================================================
# Streamlit uses session state to remember data between runs
# We initialize all our variables here so they persist

if 'reactions' not in st.session_state:
    st.session_state.reactions = []  # Stores all the reactions we've added
if 'species_list' not in st.session_state:
    st.session_state.species_list = []  # List of all unique species
if 'atomic_matrix' not in st.session_state:
    st.session_state.atomic_matrix = pd.DataFrame()  # A matrix showing atomic composition
if 'stoich_matrix' not in st.session_state:
    st.session_state.stoich_matrix = pd.DataFrame()  # Matrix showing reaction stoichiometry
if 'param_data' not in st.session_state:
    st.session_state.param_data = []  # Parameters like rate constants and pressures
if 'matrix_consistency' not in st.session_state:
    st.session_state.matrix_consistency = {'consistent': True, 'atomic_species': [], 'stoich_species': []}
if 'mass_balance_result' not in st.session_state:
    st.session_state.mass_balance_result = {'balanced': None, 'message': ''}
if 'app_mode' not in st.session_state:
    st.session_state.app_mode = 'Generator'  # Can be 'Generator' or 'Analysis'
if 'uploaded_files' not in st.session_state:
    st.session_state.uploaded_files = {'atomic': None, 'stoich': None, 'params': None}
if 'analysis_mode' not in st.session_state:
    st.session_state.analysis_mode = 'Mass Balance'  # Type of analysis in analysis mode
if 'sigma_values' not in st.session_state:
    st.session_state.sigma_values = []  # Stoichiometric numbers for thermodynamics
if 'equilibrium_constants' not in st.session_state:
    st.session_state.equilibrium_constants = []  # Calculated equilibrium constants

# =============================================================================
# HELPER FUNCTIONS FOR PARSING REACTIONS AND SPECIES
# =============================================================================

def parse_species(species_str):
    """
    Takes a species string like "2O*" or "CO(g)" and figures out 
    what the actual species is and how many of them we have.
    Returns the species name and its coefficient.
    """
    species_str = species_str.strip()  # Remove any extra spaces
    
    # Look for numbers at the beginning (like "2" in "2O*")
    coeff_match = re.match(r'^(\d+)(.+)$', species_str)
    if coeff_match:
        coeff = int(coeff_match.group(1))  # The number part
        species = coeff_match.group(2)     # The species part
    else:
        coeff = 1          # If no number, coefficient is 1
        species = species_str
    
    return species, coeff

def parse_reaction(reaction_str):
    """
    Takes a reaction string like "CO(g) + * ⇌ CO*" and breaks it down
    into reactants and products with their coefficients.
    This is the heart of the reaction parsing system!
    """
    # We only accept the equilibrium arrow ⇌ to keep things consistent
    if '⇌' not in reaction_str:
        st.error(f"❌ Please use only ⇌ arrows. Copy the arrow from the box above. Your input: {reaction_str}")
        return None, None, None
    
    # Split the reaction at the arrow
    parts = reaction_str.split('⇌')
    if len(parts) != 2:
        st.error(f"Invalid reaction format: {reaction_str}")
        return None, None, None
    
    reactants_str, products_str = parts
    
    # Parse the reactants side
    reactants = {}
    for reactant in reactants_str.split('+'):
        reactant = reactant.strip()
        if reactant:  # Skip empty strings
            species, coeff = parse_species(reactant)
            # If we see the same species twice, add up the coefficients
            reactants[species] = reactants.get(species, 0) + coeff
    
    # Parse the products side (same logic as reactants)
    products = {}
    for product in products_str.split('+'):
        product = product.strip()
        if product:
            species, coeff = parse_species(product)
            products[species] = products.get(species, 0) + coeff
    
    return reactants, products, reaction_str

# =============================================================================
# CHEMICAL FORMULA PARSING FUNCTIONS
# =============================================================================
# These functions figure out the atomic composition of any chemical species

def get_atomic_composition():
    """
    Provides a fallback dictionary for really basic species.
    Most species will be parsed automatically now, but we keep
    this for the simplest cases like the surface site "*".
    """
    compositions = {
        '*': {'*': 1},  # Empty surface site - can't parse this automatically!
    }
    return compositions

def parse_species_formula(formula):
    """
    This is our enhanced chemical formula parser! It can handle complex
    formulas like Ca(OH)2, Al2(SO4)3, and even surface species like CO*.
    Much more powerful than the old hardcoded approach.
    """
    composition = {}
    
    # Handle the special case of empty surface sites
    if formula == '*':
        return {'*': 1}
    
    # Check if this is an adsorbed species (ends with *)
    surface_site = False
    if formula.endswith('*'):
        surface_site = True
        formula = formula[:-1]  # Remove the * for now, we'll add it back later
    
    def expand_parentheses(formula_str):
        """
        Handles parentheses in chemical formulas like Ca(OH)2.
        This expands Ca(OH)2 into CaO2H2 so we can parse it easily.
        """
        import re
        
        # Keep expanding until no more parentheses
        while '(' in formula_str:
            # Find the innermost parentheses (handles nested cases)
            match = re.search(r'\(([^()]+)\)(\d*)', formula_str)
            if not match:
                break  # Something went wrong, bail out
            
            group_content = match.group(1)    # What's inside the parentheses
            multiplier = int(match.group(2)) if match.group(2) else 1  # Number after parentheses
            
            # Expand everything inside the parentheses
            expanded = ''
            group_elements = re.findall(r'([A-Z][a-z]?)(\d*)', group_content)
            for element, count in group_elements:
                count = int(count) if count else 1
                expanded += element + str(count * multiplier)
            
            # Replace the parentheses group with the expanded version
            formula_str = formula_str[:match.start()] + expanded + formula_str[match.end():]
        
        return formula_str
    
    # First, expand any parentheses we find
    expanded_formula = expand_parentheses(formula)
    
    # Now parse all the element-number pairs
    import re
    pattern = r'([A-Z][a-z]?)(\d*)'  # Matches things like "Ca", "O2", "Al", etc.
    matches = re.findall(pattern, expanded_formula)
    
    # Build up our composition dictionary
    for element, count in matches:
        count = int(count) if count else 1  # If no number, assume 1
        composition[element] = composition.get(element, 0) + count
    
    # Add the surface site back if this was an adsorbed species
    if surface_site:
        composition['*'] = 1
    
    # If we couldn't parse anything, mark it as unknown
    if not composition and formula:
        composition = {'Unknown': 1}
        print(f"Warning: Could not parse formula '{formula}'. Treating as unknown species.")
    
    return composition

def get_species_atomic_composition(species):
    """
    Main function that gets the atomic composition for any species.
    First tries the fallback dictionary, then uses our smart parser.
    """
    # Remove phase notation like (g), (l), (s) - we don't need it for composition
    clean_species = re.sub(r'\([gls]\)', '', species)
    
    # Check the basic fallback dictionary first
    compositions = get_atomic_composition()
    if clean_species in compositions:
        return compositions[clean_species]
    
    # Use our enhanced parser for everything else
    return parse_species_formula(clean_species)

# =============================================================================
# PARAMETER AND MATRIX GENERATION FUNCTIONS
# =============================================================================

def create_default_parameters():
    """
    Creates a default set of parameters for the reactions we've added.
    This includes temperature, gas constant, pressures, rate constants,
    and some placeholder constants. Users can edit these later.
    """
    params = []
    
    # Global thermodynamic parameters
    params.append({'Reaction_Descrp': '', 'Parameter': 'T', 'Values': 320.0, 'Units': 'K'})
    params.append({'Reaction_Descrp': '', 'Parameter': 'R', 'Values': 8.31446, 'Units': 'JK^-1mol^-1'})
    
    # Create pressure parameters for all gas phase species
    gas_species = [s for s in st.session_state.species_list if not s.endswith('*')]
    for i, species in enumerate(gas_species):
        params.append({'Reaction_Descrp': species, 'Parameter': f'P{i+1}', 'Values': 1.0e-8, 'Units': 'bar'})
    
    # Create forward and reverse rate constants for each reaction
    for i in range(len(st.session_state.reactions)):
        reaction_id = f'r{i+1}'
        params.append({'Reaction_Descrp': reaction_id, 'Parameter': f'k{i+1}f', 'Values': 1.0, 'Units': '-'})
        params.append({'Reaction_Descrp': '', 'Parameter': f'k{i+1}r', 'Values': 1.0, 'Units': '-'})
    
    # Add some placeholder constants (users might need these for kinetic modeling)
    for i in range(len(st.session_state.reactions)):
        for direction in ['f', 'r']:  # Forward and reverse
            for const_type in ['a', 'b', 'c']:  # Different types of constants
                params.append({'Reaction_Descrp': 'const', 'Parameter': f'{const_type}{i+1}{direction}', 'Values': 1.0, 'Units': '-'})
    
    return params

def update_matrices():
    """
    This is a big function that rebuilds both the atomic and stoichiometric
    matrices whenever we add or remove reactions. It makes sure everything
    stays consistent and properly ordered.
    """
    # If no reactions, just clear everything
    if not st.session_state.reactions:
        st.session_state.atomic_matrix = pd.DataFrame()
        st.session_state.stoich_matrix = pd.DataFrame()
        st.session_state.species_list = []
        return
    
    # Collect all unique species from all reactions
    all_species = set()
    for _, reactants, products, _ in st.session_state.reactions:
        all_species.update(reactants.keys())
        all_species.update(products.keys())
    
    # Organize species: gas phase first, then surface species, with * at the very end
    gas_species = sorted([s for s in all_species if not s.endswith('*')])
    surface_species = [s for s in all_species if s.endswith('*')]
    
    # Sort surface species but keep empty sites (*) at the end
    other_surface = sorted([s for s in surface_species if s != '*'])
    empty_sites = [s for s in surface_species if s == '*']
    surface_species = other_surface + empty_sites
    
    # Create our final ordered species list
    ordered_species = gas_species + surface_species
    st.session_state.species_list = ordered_species
    
    # ======================
    # BUILD THE ATOMIC MATRIX
    # ======================
    
    # Find all the different atoms we need to track
    all_atoms = set()
    for species in ordered_species:
        composition = get_species_atomic_composition(species)
        all_atoms.update(composition.keys())
    
    # Sort atoms with surface sites (*) at the end
    other_atoms = sorted([atom for atom in all_atoms if atom != '*'])
    surface_atom = [atom for atom in all_atoms if atom == '*']
    all_atoms = other_atoms + surface_atom
    
    # Build the atomic matrix data
    atomic_data = []
    for atom in all_atoms:
        row = [atom]  # First column is the atom name
        for species in ordered_species:
            composition = get_species_atomic_composition(species)
            row.append(composition.get(atom, 0))  # How many of this atom in this species
        atomic_data.append(row)
    
    # Create the atomic matrix DataFrame
    columns = ['A\\S'] + ordered_species  # A\S means "Atoms \ Species"
    st.session_state.atomic_matrix = pd.DataFrame(atomic_data, columns=columns)
    
    # ============================
    # BUILD THE STOICHIOMETRIC MATRIX
    # ============================
    
    stoich_data = []
    
    # Create column headers with proper notation
    gas_species = [s for s in ordered_species if not s.endswith('*')]
    surface_species = [s for s in ordered_species if s.endswith('*')]
    stoich_columns = ['r\\S'] + [f'P_{s}' for s in gas_species] + [f'theta_{s}' for s in surface_species]
    
    # Build each row (one for each reaction)
    for i, (_, reactants, products, _) in enumerate(st.session_state.reactions):
        row = [f'r{i+1}']  # First column is reaction name
        
        # Gas phase species columns
        for species in gas_species:
            # Net stoichiometric coefficient = products - reactants
            net_coeff = products.get(species, 0) - reactants.get(species, 0)
            row.append(net_coeff)
        
        # Surface species columns
        for species in surface_species:
            net_coeff = products.get(species, 0) - reactants.get(species, 0)
            row.append(net_coeff)
        
        stoich_data.append(row)
    
    # Create the stoichiometric matrix DataFrame
    st.session_state.stoich_matrix = pd.DataFrame(stoich_data, columns=stoich_columns)
    
    # Check that both matrices have consistent species ordering
    verify_matrix_consistency()

def verify_matrix_consistency():
    """
    Double-checks that our atomic and stoichiometric matrices have
    the same species in the same order. This is crucial for mass balance
    calculations to work properly.
    """
    if st.session_state.atomic_matrix.empty or st.session_state.stoich_matrix.empty:
        return True
    
    # Get species from atomic matrix (skip first column which is atom names)
    atomic_species = st.session_state.atomic_matrix.columns[1:].tolist()
    
    # Get species from stoichiometric matrix and clean up the prefixes
    stoich_species = []
    for col in st.session_state.stoich_matrix.columns[1:]:
        if col.startswith('P_'):
            stoich_species.append(col[2:])  # Remove 'P_' prefix
        elif col.startswith('theta_'):
            stoich_species.append(col[6:])  # Remove 'theta_' prefix
    
    # Check if the lists match
    consistency_check = atomic_species == stoich_species
    
    # Store the result for display purposes
    st.session_state.matrix_consistency = {
        'consistent': consistency_check,
        'atomic_species': atomic_species,
        'stoich_species': stoich_species
    }
    
    return consistency_check

# =============================================================================
# ANALYSIS AND CALCULATION FUNCTIONS
# =============================================================================

def check_mass_balance(atomic_matrix=None, stoich_matrix=None):
    """
    Performs mass balance checking using the fundamental equation A*v=0,
    where A is the atomic matrix and v is the stoichiometric vector.
    This checks if atoms are conserved across all reactions.
    """
    # Use provided matrices or default to session state
    if atomic_matrix is None:
        atomic_matrix = st.session_state.atomic_matrix
    if stoich_matrix is None:
        stoich_matrix = st.session_state.stoich_matrix
        
    if atomic_matrix.empty or stoich_matrix.empty:
        return False, "Matrices not available for mass balance check"
    
    # Extract the numerical parts (skip label columns)
    atomic_matrix_vals = atomic_matrix.iloc[:, 1:].values
    stoich_matrix_vals = stoich_matrix.iloc[:, 1:].values
    
    errors = []
    mass_balanced = True
    
    # Check each reaction individually
    for i, reaction_row in enumerate(stoich_matrix_vals):
        # Calculate A * v for this reaction
        mass_balance_result = np.dot(atomic_matrix_vals, reaction_row)
        
        # If any element is not close to zero, we have an imbalance
        if not np.allclose(mass_balance_result, 0, atol=1e-10):
            mass_balanced = False
            reaction_name = stoich_matrix.iloc[i, 0]
            errors.append({
                'reaction': i + 1,
                'reaction_name': reaction_name,
                'imbalance': mass_balance_result.tolist()
            })
    
    # Prepare the result message
    if mass_balanced:
        return True, "✅ Mass is conserved across all reactions."
    else:
        error_message = "❌ Mass is NOT conserved in the following reactions:\n"
        for error in errors:
            error_message += f"• Reaction {error['reaction']} ({error['reaction_name']}): Imbalance = {error['imbalance']}\n"
        error_message += "\nPlease check and correct the Atomic and/or Stoichiometric matrices."
        return False, error_message

def calculate_equilibrium_constants(params_df, sigma_values, k_pattern="auto"):
    """
    Calculates individual and overall equilibrium constants from rate constants.
    The overall equilibrium constant is calculated as K_overall = ∏K_i^σ_i
    where K_i = k_forward/k_reverse for each reaction.
    """
    try:
        # Extract temperature and gas constant from parameters
        T_row = params_df[params_df['Parameter'] == 'T']
        R_row = params_df[params_df['Parameter'] == 'R']
        
        # Use defaults if not found
        if T_row.empty:
            st.warning("Temperature (T) not found, using default 320 K")
            T = 320.0
        else:
            T = T_row['Values'].iloc[0]
            
        if R_row.empty:
            st.warning("Gas constant (R) not found, using default 8.314 J/(mol·K)")
            R = 8.314
        else:
            R = R_row['Values'].iloc[0]
        
        # Find all rate constants - we need to be flexible with naming conventions
        k_forward = []
        k_reverse = []
        
        # Get all parameters that look like rate constants
        k_params = params_df[params_df['Parameter'].str.contains(r'^k\d+', regex=True, na=False)]
        k_params = k_params.sort_values('Parameter')
        
        if k_pattern == "alternating":
            # User told us every other k is reverse (k1=forward, k2=reverse, etc.)
            k_values = k_params['Values'].tolist()
            for i in range(0, len(k_values), 2):
                if i < len(k_values):
                    k_forward.append(k_values[i])
                if i + 1 < len(k_values):
                    k_reverse.append(k_values[i + 1])
                else:
                    k_reverse.append(1.0)  # Default if missing
                    
        else:
            # Try to automatically detect the naming pattern
            for i in range(len(sigma_values)):
                # Try different possible naming conventions
                possible_forward = [f'k{i+1}f', f'k{i+1}F', f'k_{i+1}_f', f'k{i+1}', f'k{2*i+1}']
                possible_reverse = [f'k{i+1}r', f'k{i+1}R', f'k_{i+1}_r', f'k{i+1}_rev', f'k{2*i+2}']
                
                kf_found = False
                kr_found = False
                
                # Look for forward rate constant
                for kf_name in possible_forward:
                    kf_row = params_df[params_df['Parameter'] == kf_name]
                    if not kf_row.empty:
                        k_forward.append(kf_row['Values'].iloc[0])
                        kf_found = True
                        break
                
                # Look for reverse rate constant
                for kr_name in possible_reverse:
                    kr_row = params_df[params_df['Parameter'] == kr_name]
                    if not kr_row.empty:
                        k_reverse.append(kr_row['Values'].iloc[0])
                        kr_found = True
                        break
                
                # Use defaults if we couldn't find them
                if not kf_found:
                    st.warning(f"Forward rate constant for reaction {i+1} not found, using 1.0")
                    k_forward.append(1.0)
                    
                if not kr_found:
                    st.warning(f"Reverse rate constant for reaction {i+1} not found, using 1.0")
                    k_reverse.append(1.0)
        
        # Make sure we have the right number of rate constants
        while len(k_forward) < len(sigma_values):
            k_forward.append(1.0)
        while len(k_reverse) < len(sigma_values):
            k_reverse.append(1.0)
            
        # Calculate individual equilibrium constants
        K_eq = []
        for i in range(len(sigma_values)):
            if k_reverse[i] != 0:
                K_i = k_forward[i] / k_reverse[i]
                K_eq.append(K_i)
            else:
                st.warning(f"Zero reverse rate constant for reaction {i+1}")
                K_eq.append(1.0)
        
        # Calculate overall equilibrium constant: K_overall = ∏(K_i^σ_i)
        K_overall = 1.0
        for i, (K_i, sigma_i) in enumerate(zip(K_eq, sigma_values)):
            if K_i > 0:  # Make sure we can take the logarithm if needed
                K_overall *= K_i ** sigma_i
            else:
                st.warning(f"Non-positive equilibrium constant for reaction {i+1}")
                K_overall *= 1.0
        
        return K_eq, K_overall, T, R, k_forward, k_reverse
        
    except Exception as e:
        st.error(f"Error calculating equilibrium constants: {e}")
        return [], 1.0, 300.0, 8.314, [], []

# =============================================================================
# FILE HANDLING FUNCTIONS
# =============================================================================

def load_uploaded_files():
    """
    Loads and validates CSV files uploaded by the user in Analysis mode.
    Returns a dictionary with the loaded DataFrames or None if there were errors.
    """
    files_loaded = {}
    
    # Try to load atomic matrix
    if st.session_state.uploaded_files['atomic'] is not None:
        try:
            atomic_df = pd.read_csv(st.session_state.uploaded_files['atomic'])
            files_loaded['atomic'] = atomic_df
        except Exception as e:
            st.error(f"Error loading atomic matrix: {e}")
            return None
    
    # Try to load stoichiometric matrix
    if st.session_state.uploaded_files['stoich'] is not None:
        try:
            stoich_df = pd.read_csv(st.session_state.uploaded_files['stoich'])
            files_loaded['stoich'] = stoich_df
        except Exception as e:
            st.error(f"Error loading stoichiometric matrix: {e}")
            return None
    
    # Try to load parameters
    if st.session_state.uploaded_files['params'] is not None:
        try:
            params_df = pd.read_csv(st.session_state.uploaded_files['params'])
            files_loaded['params'] = params_df
        except Exception as e:
            st.error(f"Error loading parameters: {e}")
            return None
    
    return files_loaded

# =============================================================================
# MAIN APPLICATION INTERFACE
# =============================================================================

# Mode selector in the sidebar - lets user choose between Generator and Analysis
with st.sidebar:
    st.header("Application Mode")
    mode = st.radio(
        "Select Mode:",
        ["Generator", "Analysis"],
        index=0 if st.session_state.app_mode == 'Generator' else 1
    )
    st.session_state.app_mode = mode

# =============================================================================
# GENERATOR MODE - BUILD MECHANISMS FROM SCRATCH
# =============================================================================

if st.session_state.app_mode == 'Generator':
    # Sidebar controls for adding reactions
    with st.sidebar:
        st.header("Reaction Input")
        
        # Reminder about proper formatting
        st.markdown("**Use ⇌ arrow (copy above) • One reaction at a time**")
        
        reaction_input = st.text_area(
            "Enter ONE reaction:",
            placeholder="CO(g) + * ⇌ CO*",
            help="Use ⇌ for arrows. Use * for surface sites. Add ONE reaction, then click 'Add Reaction'.",
            height=80
        )
        
        # Button to add the reaction
        if st.button("Add Reaction", type="primary"):
            if reaction_input.strip():
                reactants, products, original = parse_reaction(reaction_input.strip())
                if reactants is not None and products is not None:
                    # Add the reaction to our list
                    st.session_state.reactions.append((len(st.session_state.reactions), reactants, products, original))
                    # Rebuild matrices and parameters
                    update_matrices()
                    st.session_state.param_data = create_default_parameters()
                    st.session_state.mass_balance_result = {'balanced': None, 'message': ''}
                    st.success("Reaction added!")
                    st.rerun()
        
        # Show current reactions with delete buttons
        if st.session_state.reactions:
            st.subheader("Current Reactions")
            for i, (idx, reactants, products, original) in enumerate(st.session_state.reactions):
                col1, col2 = st.columns([3, 1])
                with col1:
                    st.write(f"r{i+1}: {original}")
                with col2:
                    if st.button("❌", key=f"del_{i}", help="Delete reaction"):
                        st.session_state.reactions.pop(i)
                        update_matrices()
                        st.session_state.param_data = create_default_parameters()
                        st.session_state.mass_balance_result = {'balanced': None, 'message': ''}
                        st.rerun()
        
        # Clear all button
        if st.button("Clear All", type="secondary"):
            st.session_state.reactions = []
            st.session_state.species_list = []
            st.session_state.atomic_matrix = pd.DataFrame()
            st.session_state.stoich_matrix = pd.DataFrame()
            st.session_state.param_data = []
            st.session_state.mass_balance_result = {'balanced': None, 'message': ''}
            st.rerun()

    # Main content area - only shows if we have reactions
    if st.session_state.reactions:
        # Create tabs for different views
        tab1, tab2, tab3, tab4, tab5 = st.tabs(["📋 Reaction Mechanism", "📊 Stoichiometric Matrix", "⚛️ Atomic Matrix", "⚙️ Parameters", "🌡️ Thermodynamics"])
        
        # =============================
        # TAB 1: REACTION MECHANISM VIEW
        # =============================
        with tab1:
            st.subheader("Reaction Mechanism")
            
            # Display each reaction step nicely
            st.markdown("**Current Reaction Steps:**")
            
            for i, (_, reactants, products, original) in enumerate(st.session_state.reactions):
                # Create a clean layout for each reaction
                col1, col2 = st.columns([1, 8])
                
                with col1:
                    st.markdown(f"**Step {i+1}:**")
                
                with col2:
                    st.code(f"r{i+1}: {original}", language=None)
            
            # Show summary info if we have multiple reactions
            if len(st.session_state.reactions) > 1:
                st.markdown("---")
                st.markdown(f"**Total Steps:** {len(st.session_state.reactions)}")
        
        # ===================================
        # TAB 2: STOICHIOMETRIC MATRIX
        # ===================================
        with tab2:
            st.subheader("Stoichiometric Matrix")
            if not st.session_state.stoich_matrix.empty:
                st.dataframe(st.session_state.stoich_matrix, use_container_width=True)
                
                # Download section
                col1, col2 = st.columns([2, 1])
                with col1:
                    stoich_filename = st.text_input("Stoichiometric matrix filename:", value="Stoich_1.csv")
                with col2:
                    csv_stoich = st.session_state.stoich_matrix.to_csv(index=False)
                    st.download_button(
                        label="📥 Download CSV",
                        data=csv_stoich,
                        file_name=stoich_filename,
                        mime="text/csv"
                    )
        
        # ==========================
        # TAB 3: ATOMIC MATRIX
        # ==========================
        with tab3:
            st.subheader("Atomic Matrix")
            if not st.session_state.atomic_matrix.empty:
                st.dataframe(st.session_state.atomic_matrix, use_container_width=True)
                
                # Download section
                col1, col2 = st.columns([2, 1])
                with col1:
                    atomic_filename = st.text_input("Atomic matrix filename:", value="Atomic_1.csv")
                with col2:
                    csv_atomic = st.session_state.atomic_matrix.to_csv(index=False)
                    st.download_button(
                        label="📥 Download CSV",
                        data=csv_atomic,
                        file_name=atomic_filename,
                        mime="text/csv"
                    )
        
        # =======================
        # TAB 4: PARAMETERS
        # =======================
        with tab4:
            st.subheader("Parameter Table")
            
            if st.session_state.param_data:
                # Create an editable table of parameters
                param_df = pd.DataFrame(st.session_state.param_data)
                
                # Display the editable table
                edited_df = st.data_editor(
                    param_df,
                    use_container_width=True,
                    num_rows="dynamic",  # Allows adding/removing rows
                    column_config={
                        "Reaction_Descrp": st.column_config.TextColumn("Reaction Description"),
                        "Parameter": st.column_config.TextColumn("Parameter"),
                        "Values": st.column_config.NumberColumn("Values", format="%.2e"),
                        "Units": st.column_config.TextColumn("Units")
                    }
                )
                
                # Update our session state with any edits
                st.session_state.param_data = edited_df.to_dict('records')
                
                # Download section
                col1, col2 = st.columns([2, 1])
                with col1:
                    param_filename = st.text_input("Parameter file filename:", value="Param_1.csv")
                with col2:
                    csv_param = edited_df.to_csv(index=False)
                    st.download_button(
                        label="📥 Download CSV",
                        data=csv_param,
                        file_name=param_filename,
                        mime="text/csv"
                    )
        
        # ==============================
        # TAB 5: THERMODYNAMICS ANALYSIS
        # ==============================
        with tab5:
            st.subheader("Thermodynamics Analysis")
            
            if st.session_state.param_data:
                # Convert our parameter data to a DataFrame for processing
                params_df = pd.DataFrame(st.session_state.param_data)
                
                # Parameter inspection section - helps users understand their data
                with st.expander("🔍 Inspect Generated Parameters"):
                    st.write("**Current Parameters:**")
                    st.dataframe(params_df, use_container_width=True)
                    
                    # Highlight rate constants specifically
                    k_params = params_df[params_df['Parameter'].str.contains(r'^k\d+', regex=True, na=False)]
                    if not k_params.empty:
                        st.write("**Rate Constants Found:**")
                        st.dataframe(k_params, use_container_width=True)
                
                # Let user choose how their rate constants are organized
                st.markdown("**Rate Constant Pattern:**")
                k_pattern_gen = st.radio(
                    "How are your rate constants organized?",
                    ["auto", "alternating"],
                    format_func=lambda x: {
                        "auto": "Auto-detect (kXf/kXr pattern)",
                        "alternating": "Alternating (k1=forward, k2=reverse, k3=forward, k4=reverse...)"
                    }[x],
                    help="Select how your rate constants are named/organized",
                    key="gen_k_pattern"
                )
                
                # Allow inline parameter editing
                st.markdown("**Edit Parameters:**")
                edited_params_gen = st.data_editor(
                    params_df,
                    use_container_width=True,
                    num_rows="dynamic",
                    column_config={
                        "Reaction_Descrp": st.column_config.TextColumn("Reaction Description"),
                        "Parameter": st.column_config.TextColumn("Parameter"),
                        "Values": st.column_config.NumberColumn("Values", format="%.6e"),
                        "Units": st.column_config.TextColumn("Units")
                    },
                    key="gen_params_editor"
                )
                
                # Update session state with any parameter edits
                st.session_state.param_data = edited_params_gen.to_dict('records')
                
                # Get the number of reactions for sigma input
                num_reactions = len(st.session_state.reactions)
                
                st.markdown("---")
                st.subheader("Stoichiometric Numbers (σ)")
                st.markdown("Enter the stoichiometric number for each reaction step:")
                
                # Initialize sigma values for generator mode
                if f'sigma_values_gen' not in st.session_state:
                    st.session_state.sigma_values_gen = [1.0] * num_reactions
                elif len(st.session_state.sigma_values_gen) != num_reactions:
                    st.session_state.sigma_values_gen = [1.0] * num_reactions
                
                # Create input fields for sigma values
                sigma_inputs_gen = []
                cols = st.columns(min(4, num_reactions))
                
                for i in range(num_reactions):
                    col_idx = i % 4
                    with cols[col_idx]:
                        sigma_val = st.number_input(
                            f"σ{i+1}",
                            value=st.session_state.sigma_values_gen[i],
                            step=0.1,
                            format="%.3f",
                            key=f"gen_sigma_{i}",
                            help=f"Stoichiometric number for reaction r{i+1}"
                        )
                        sigma_inputs_gen.append(sigma_val)
                
                st.session_state.sigma_values_gen = sigma_inputs_gen
                
                # Show current reactions for reference
                with st.expander("📋 Current Reactions"):
                    for i, (_, reactants, products, original) in enumerate(st.session_state.reactions):
                        st.write(f"**r{i+1}:** {original}")
                
                # Main calculation button
                if st.button("🧮 Calculate Overall Equilibrium", type="primary", key="gen_calc_eq"):
                    K_eq_list_gen, K_overall_gen, T_gen, R_gen, k_forward_gen, k_reverse_gen = calculate_equilibrium_constants(
                        edited_params_gen, st.session_state.sigma_values_gen, k_pattern_gen
                    )
                    
                    st.markdown("---")
                    st.subheader("🎯 Results")
                    
                    # Show which rate constants were actually used
                    st.markdown("**Rate Constants Used:**")
                    rate_data_gen = []
                    for i, (kf, kr) in enumerate(zip(k_forward_gen, k_reverse_gen)):
                        # Get the original reaction for reference
                        original_reaction = st.session_state.reactions[i][3] if i < len(st.session_state.reactions) else f"r{i+1}"
                        rate_data_gen.append({
                            'Reaction': f'r{i+1}',
                            'Original': original_reaction,
                            'k_forward': f'{kf:.6e}',
                            'k_reverse': f'{kr:.6e}',
                            'K_eq (kf/kr)': f'{kf/kr if kr != 0 else "∞":.6e}'
                        })
                    
                    rate_df_gen = pd.DataFrame(rate_data_gen)
                    st.dataframe(rate_df_gen, use_container_width=True)
                    
                    # Show the equilibrium calculation breakdown
                    st.markdown("**Equilibrium Calculation:**")
                    eq_data_gen = []
                    for i, (K_i, sigma_i) in enumerate(zip(K_eq_list_gen, st.session_state.sigma_values_gen)):
                        eq_data_gen.append({
                            'Reaction': f'r{i+1}',
                            'K_eq': f'{K_i:.6e}',
                            'σ': f'{sigma_i:.3f}',
                            'K^σ': f'{K_i**sigma_i:.6e}'
                        })
                    
                    eq_df_gen = pd.DataFrame(eq_data_gen)
                    st.dataframe(eq_df_gen, use_container_width=True)
                    
                    # Display the overall equilibrium constant
                    st.markdown("**Overall Equilibrium Constant:**")
                    st.latex(r"K_{overall,eq} = \prod_{i} K_{i,eq}^{\sigma_i}")
                    
                    # Format the result nicely in scientific notation
                    if K_overall_gen != 0:
                        exponent = math.floor(math.log10(abs(K_overall_gen)))
                        mantissa = K_overall_gen / (10 ** exponent)
                        st.markdown(f"**K_overall = {mantissa:.4f} × 10^{exponent}**")
                    else:
                        st.markdown("**K_overall = 0**")
                    
                    # Additional thermodynamic information
                    col1, col2, col3 = st.columns(3)
                    with col1:
                        st.metric("Temperature", f"{T_gen:.1f} K")
                    with col2:
                        st.metric("Gas Constant", f"{R_gen:.5f} J/(mol·K)")
                    with col3:
                        if K_overall_gen > 0:
                            ln_K = math.log(K_overall_gen)
                            st.metric("ln(K_overall)", f"{ln_K:.4f}")
                
            else:
                st.info("Add reactions first to generate parameters for thermodynamics analysis.")
    
        # =============================
        # CONSISTENCY AND MASS BALANCE CHECKS
        # =============================
        
        st.markdown("---")
        
        # Matrix consistency check
        if 'matrix_consistency' in st.session_state:
            if st.session_state.matrix_consistency['consistent']:
                st.success("✅ Matrix Consistency Check: Both atomic and stoichiometric matrices have consistent species ordering (Gas → Surface)")
            else:
                st.error("❌ Matrix Consistency Check: Species ordering mismatch detected!")
                with st.expander("View Ordering Details"):
                    col1, col2 = st.columns(2)
                    with col1:
                        st.write("**Atomic Matrix Species Order:**")
                        for i, species in enumerate(st.session_state.matrix_consistency['atomic_species']):
                            st.write(f"{i+1}. {species}")
                    with col2:
                        st.write("**Stoichiometric Matrix Species Order:**")
                        for i, species in enumerate(st.session_state.matrix_consistency['stoich_species']):
                            st.write(f"{i+1}. {species}")
        
        # Mass balance verification section
        st.markdown("### 🧮 Mass Balance Verification")
        col1, col2 = st.columns([1, 3])
        
        with col1:
            if st.button("Check Mass Balance", type="primary"):
                balanced, message = check_mass_balance()
                st.session_state.mass_balance_result = {'balanced': balanced, 'message': message}
                st.rerun()
        
        with col2:
            if st.session_state.mass_balance_result['balanced'] is not None:
                if st.session_state.mass_balance_result['balanced']:
                    st.success(st.session_state.mass_balance_result['message'])
                else:
                    st.error(st.session_state.mass_balance_result['message'])
        
        # Information about mass balance
        with st.expander("ℹ️ About Mass Balance Check"):
            st.markdown("""
            **Mass Balance Verification** ensures that atoms are conserved across all reactions.
            
            **How it works:**
            - For each reaction, the atomic matrix (A) is multiplied by the stoichiometric coefficients (v)
            - If A × v = 0 for all reactions, mass is conserved
            - Any non-zero result indicates a mass imbalance
            
            **Example:**
            For reaction: CO + O* → CO2 + *
            - Carbon: 1×1 + 0×1 = 1 (reactants) vs 1×1 + 0×1 = 1 (products) ✓
            - Oxygen: 1×1 + 1×1 = 2 (reactants) vs 2×1 + 0×1 = 2 (products) ✓
            """)
        
        # Summary metrics
        col1, col2, col3, col4 = st.columns(4)
        with col1:
            st.metric("Reactions", len(st.session_state.reactions))
        with col2:
            st.metric("Species", len(st.session_state.species_list))
        with col3:
            gas_count = len([s for s in st.session_state.species_list if not s.endswith('*')])
            st.metric("Gas Species", gas_count)
        with col4:
            surface_count = len([s for s in st.session_state.species_list if s.endswith('*')])
            st.metric("Surface Species", surface_count)
        
        # Species ordering display
        with st.expander("📋 Current Species Ordering (Gas → Surface, * Last)"):
            col1, col2 = st.columns(2)
            gas_species = [s for s in st.session_state.species_list if not s.endswith('*')]
            surface_species = [s for s in st.session_state.species_list if s.endswith('*')]
            
            with col1:
                st.write("**Gas Species:**")
                for i, species in enumerate(gas_species):
                    st.write(f"{i+1}. {species}")
            
            with col2:
                st.write("**Surface Species (Empty Sites Last):**")
                for i, species in enumerate(surface_species):
                    if species == '*':
                        st.write(f"{i+1}. {species} ← Empty Sites")
                    else:
                        st.write(f"{i+1}. {species}")

    else:
        # Welcome screen when no reactions have been added yet
        st.markdown("""
        ## 🚀 Getting Started
        
        Enter chemical reactions **one at a time** in the sidebar to automatically generate:
        - **Stoichiometric Matrix**: Shows how each reaction affects each species
        - **Atomic Matrix**: Shows atomic composition of each species  
        - **Parameter Table**: Editable parameters for reaction kinetics
        - **Thermodynamics Analysis**: Calculate overall equilibrium constants
        
        All matrices maintain consistent ordering: **Gas → Surface → Empty Sites (*)** last.
        
        ### ⚠️ Important Rules:
        - **Use ONLY ⇌ arrows** (copy from the box above)
        - **Enter ONE reaction at a time**, then click "Add Reaction"
        - **Use `*` for surface sites**
        - **Use `(g)` for gas phase** (optional)
        
        ### 📝 Example Reactions (enter one by one):
        ```
        CO(g) + * ⇌ CO*
        O2(g) + * ⇌ O2*
        O2* + * ⇌ 2O*
        CO* + O* ⇌ CO2(g) + 2*
        ```
        
        ### 💡 Features:
        - **Automatic Ordering**: Gas species columns appear first, then surface species
        - **Empty Sites Last**: `*` (empty sites) always appear in the last column/row
        - **Consistency Check**: Both matrices maintain identical species ordering
        - **Mass Balance**: Verify atom conservation across all reactions
        - **Thermodynamics**: Calculate K_overall = ∏K_i^σ_i using generated parameters
        
        ### 🌡️ Thermodynamics Analysis:
        After adding reactions, use the **Thermodynamics** tab to:
        - Edit rate constants and parameters
        - Set stoichiometric numbers (σ) for each step
        - Calculate individual and overall equilibrium constants
        """)

# =============================================================================
# ANALYSIS MODE - ANALYZE UPLOADED FILES
# =============================================================================

else:
    # Sidebar for file uploads
    with st.sidebar:
        st.header("File Upload")
        
        # File upload widgets
        atomic_file = st.file_uploader(
            "Upload Atomic Matrix CSV",
            type=['csv'],
            key="atomic_upload"
        )
        if atomic_file is not None:
            st.session_state.uploaded_files['atomic'] = atomic_file
        
        stoich_file = st.file_uploader(
            "Upload Stoichiometric Matrix CSV",
            type=['csv'],
            key="stoich_upload"
        )
        if stoich_file is not None:
            st.session_state.uploaded_files['stoich'] = stoich_file
        
        params_file = st.file_uploader(
            "Upload Parameters CSV",
            type=['csv'],
            key="params_upload"
        )
        if params_file is not None:
            st.session_state.uploaded_files['params'] = params_file
        
        st.markdown("---")
        
        # Analysis type selector
        st.header("Analysis Type")
        analysis_mode = st.radio(
            "Select Analysis:",
            ["Mass Balance", "Thermodynamics"],
            index=0 if st.session_state.analysis_mode == 'Mass Balance' else 1
        )
        st.session_state.analysis_mode = analysis_mode

    # Main content for Analysis mode
    st.markdown("## 🔬 Analysis Mode")
    
    # Check if any files have been uploaded
    files_uploaded = any(st.session_state.uploaded_files[key] is not None for key in st.session_state.uploaded_files)
    
    if not files_uploaded:
        # Show instructions if no files uploaded
        st.info("Please upload your files in the sidebar to begin analysis.")
        st.markdown("""
        ### 📁 Required Files:
        1. **Atomic Matrix CSV**: Contains atomic composition of species
        2. **Stoichiometric Matrix CSV**: Contains reaction coefficients
        3. **Parameters CSV**: Contains rate constants and other parameters
        
        ### 🔍 Available Analysis Types:
        - **Mass Balance**: Verify atom conservation across reactions
        - **Thermodynamics**: Calculate overall equilibrium constant using K_overall = ∏K_i^σ_i
        """)
    else:
        # Load the uploaded files
        loaded_files = load_uploaded_files()
        
        if loaded_files:
            # Show success and file info
            st.success("✅ Files loaded successfully!")
            
            col1, col2, col3 = st.columns(3)
            with col1:
                if 'atomic' in loaded_files:
                    st.metric("Atomic Matrix", f"{loaded_files['atomic'].shape[0]}×{loaded_files['atomic'].shape[1]}")
            with col2:
                if 'stoich' in loaded_files:
                    st.metric("Stoichiometric Matrix", f"{loaded_files['stoich'].shape[0]}×{loaded_files['stoich'].shape[1]}")
            with col3:
                if 'params' in loaded_files:
                    st.metric("Parameters", f"{loaded_files['params'].shape[0]} entries")
            
            # Mass Balance Analysis
            if st.session_state.analysis_mode == 'Mass Balance':
                st.markdown("### 🧮 Mass Balance Analysis")
                
                if 'atomic' in loaded_files and 'stoich' in loaded_files:
                    # Display both matrices side by side
                    col1, col2 = st.columns(2)
                    
                    with col1:
                        st.subheader("Atomic Matrix")
                        st.dataframe(loaded_files['atomic'], use_container_width=True)
                    
                    with col2:
                        st.subheader("Stoichiometric Matrix")
                        st.dataframe(loaded_files['stoich'], use_container_width=True)
                    
                    # Mass balance check button
                    if st.button("🔍 Check Mass Balance", type="primary"):
                        balanced, message = check_mass_balance(loaded_files['atomic'], loaded_files['stoich'])
                        
                        if balanced:
                            st.success(message)
                        else:
                            st.error(message)
                
                else:
                    st.warning("Both Atomic and Stoichiometric matrices are required for mass balance analysis.")
            
            # Thermodynamics Analysis
            else:
                st.markdown("### 🌡️ Thermodynamics Analysis")
                
                if 'params' in loaded_files:
                    # Display and allow editing of parameters
                    st.subheader("Parameters")
                    
                    # Parameter inspection section
                    with st.expander("🔍 Inspect Parameter File"):
                        st.write("**All Parameters in File:**")
                        param_preview = loaded_files['params'].copy()
                        
                        # Highlight rate constants
                        k_params = param_preview[param_preview['Parameter'].str.contains(r'^k\d+', regex=True, na=False)]
                        if not k_params.empty:
                            st.write("**Rate Constants Found:**")
                            st.dataframe(k_params, use_container_width=True)
                            
                            st.write("**Rate Constant Pattern Detection:**")
                            k_names = k_params['Parameter'].tolist()
                            
                            # Try to detect the naming pattern
                            has_f_suffix = any('f' in k.lower() for k in k_names)
                            has_r_suffix = any('r' in k.lower() for k in k_names)
                            has_alternating = len(k_names) >= 2
                            
                            if has_f_suffix and has_r_suffix:
                                st.success("✅ Detected: Forward/Reverse suffixes (kXf, kXr)")
                            elif has_alternating:
                                st.info("ℹ️ Detected: Sequential numbering (possibly alternating)")
                            else:
                                st.warning("⚠️ Unknown rate constant pattern")
                    
                    # Rate constant pattern selector
                    st.markdown("**Rate Constant Pattern:**")
                    k_pattern = st.radio(
                        "How are your rate constants organized?",
                        ["auto", "alternating"],
                        format_func=lambda x: {
                            "auto": "Auto-detect (kXf/kXr or similar)",
                            "alternating": "Alternating (k1=forward, k2=reverse, k3=forward, k4=reverse...)"
                        }[x],
                        help="Select how your rate constants are named/organized"
                    )
                    
                    # Editable parameters table
                    edited_params = st.data_editor(
                        loaded_files['params'],
                        use_container_width=True,
                        num_rows="dynamic",
                        column_config={
                            "Values": st.column_config.NumberColumn("Values", format="%.6e")
                        }
                    )
                    
                    # Determine number of reactions from stoichiometric matrix
                    if 'stoich' in loaded_files:
                        num_reactions = loaded_files['stoich'].shape[0]
                        
                        st.markdown("---")
                        st.subheader("Stoichiometric Numbers (σ)")
                        st.markdown("Enter the stoichiometric number for each reaction step:")
                        
                        # Initialize sigma values if needed
                        if len(st.session_state.sigma_values) != num_reactions:
                            st.session_state.sigma_values = [1.0] * num_reactions
                        
                        # Create input fields for sigma values
                        sigma_inputs = []
                        cols = st.columns(min(4, num_reactions))
                        
                        for i in range(num_reactions):
                            col_idx = i % 4
                            with cols[col_idx]:
                                sigma_val = st.number_input(
                                    f"σ{i+1}",
                                    value=st.session_state.sigma_values[i],
                                    step=0.1,
                                    format="%.3f",
                                    key=f"sigma_{i}"
                                )
                                sigma_inputs.append(sigma_val)
                        
                        st.session_state.sigma_values = sigma_inputs
                        
                        # Calculate equilibrium constants
                        if st.button("🧮 Calculate Overall Equilibrium", type="primary"):
                            K_eq_list, K_overall, T, R, k_forward, k_reverse = calculate_equilibrium_constants(
                                edited_params, st.session_state.sigma_values, k_pattern
                            )
                            st.session_state.equilibrium_constants = K_eq_list
                            
                            st.markdown("---")
                            st.subheader("🎯 Results")
                            
                            # Display rate constants used
                            st.markdown("**Rate Constants Used:**")
                            rate_data = []
                            for i, (kf, kr) in enumerate(zip(k_forward, k_reverse)):
                                rate_data.append({
                                    'Reaction': f'r{i+1}',
                                    'k_forward': f'{kf:.6e}',
                                    'k_reverse': f'{kr:.6e}',
                                    'K_eq (kf/kr)': f'{kf/kr if kr != 0 else "∞":.6e}'
                                })
                            
                            rate_df = pd.DataFrame(rate_data)
                            st.dataframe(rate_df, use_container_width=True)
                            
                            # Display equilibrium calculation
                            st.markdown("**Equilibrium Calculation:**")
                            eq_data = []
                            for i, (K_i, sigma_i) in enumerate(zip(K_eq_list, st.session_state.sigma_values)):
                                eq_data.append({
                                    'Reaction': f'r{i+1}',
                                    'K_eq': f'{K_i:.6e}',
                                    'σ': f'{sigma_i:.3f}',
                                    'K^σ': f'{K_i**sigma_i:.6e}'
                                })
                            
                            eq_df = pd.DataFrame(eq_data)
                            st.dataframe(eq_df, use_container_width=True)
                            
                            # Display overall equilibrium constant
                            st.markdown("**Overall Equilibrium Constant:**")
                            st.latex(r"K_{overall,eq} = \prod_{i} K_{i,eq}^{\sigma_i}")
                            
                            # Format scientific notation nicely
                            if K_overall != 0:
                                exponent = math.floor(math.log10(abs(K_overall)))
                                mantissa = K_overall / (10 ** exponent)
                                st.markdown(f"**K_overall = {mantissa:.4f} × 10^{exponent}**")
                            else:
                                st.markdown("**K_overall = 0**")
                            
                            # Additional thermodynamic info
                            col1, col2, col3 = st.columns(3)
                            with col1:
                                st.metric("Temperature", f"{T:.1f} K")
                            with col2:
                                st.metric("Gas Constant", f"{R:.5f} J/(mol·K)")
                            with col3:
                                if K_overall > 0:
                                    ln_K = math.log(K_overall)
                                    st.metric("ln(K_overall)", f"{ln_K:.4f}")
                    
                    else:
                        st.warning("Stoichiometric matrix is required to determine the number of reactions.")
                
                else:
                    st.warning("Parameters file is required for thermodynamics analysis.")
