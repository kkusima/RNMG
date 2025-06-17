# Author: Kenneth Kusima
# Version 1.0
# 06/2025
# =============================================================================
# REACTION NETWORK MATRIX GENERATOR
# =============================================================================
# This app helps create stoichiometric and atomic matrices from chemical 
# reaction mechanisms

# Import all the libraries we need
import streamlit as st
import pandas as pd
import numpy as np
import re
from io import StringIO
from collections import defaultdict, Counter
import math
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots

# =============================================================================
# PAGE SETUP AND CONFIGURATION
# =============================================================================

# Set up the page layout and basic configuration
st.set_page_config(
    page_title="RNMG - Reaction Network Matrix Generator",
    page_icon="🔗",
    layout="wide"  # Uses full width of browser
)

# Main title and description
st.title("RNMG - Reaction Network Matrix Generator")
st.markdown("Create stoichiometric and atomic matrices from chemical reaction mechanisms with comprehensive thermodynamic analysis")

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
if 'equilibrium_results' not in st.session_state:
    st.session_state.equilibrium_results = {}  # Equilibrium analysis results

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
    This is our comprehensive chemical formula parser! It can handle complex
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
    
    # Use our comprehensive parser for everything else
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
    Comprehensive function to calculate equilibrium constants from rate constants with comprehensive error handling.
    This function calculates individual K_eq values and the overall equilibrium constant using K_overall = ∏K_i^σ_i
    
    Args:
        params_df: DataFrame with parameters including rate constants
        sigma_values: List of stoichiometric numbers for each reaction
        k_pattern: How rate constants are organized ("auto" or "alternating")
    
    Returns:
        Tuple of (K_eq_list, K_overall, temperature, gas_constant, k_forward, k_reverse)
    """
    try:
        # Extract temperature and gas constant from parameters with better error handling
        T_row = params_df[params_df['Parameter'].str.lower() == 't']
        R_row = params_df[params_df['Parameter'].str.lower() == 'r']
        
        # Use defaults if not found, with user warnings
        if T_row.empty:
            st.warning("Temperature (T) not found in parameters, using default 320 K")
            T = 320.0
        else:
            T = float(T_row['Values'].iloc[0])
            
        if R_row.empty:
            st.warning("Gas constant (R) not found in parameters, using default 8.314 J/(mol·K)")
            R = 8.314
        else:
            R = float(R_row['Values'].iloc[0])
        
        # Find all rate constants with improved pattern matching
        k_forward = []
        k_reverse = []
        
        # Get all parameters that look like rate constants
        k_params = params_df[params_df['Parameter'].str.contains(r'k', regex=True, na=False, case=False)]
        k_params = k_params.sort_values('Parameter')
        
        if k_pattern == "alternating":
            # Every other k is reverse (k1=forward, k2=reverse, etc.)
            k_values = []
            for _, row in k_params.iterrows():
                param_name = str(row['Parameter'])
                if 'k' in param_name.lower():
                    try:
                        k_values.append(float(row['Values']))
                    except (ValueError, TypeError):
                        st.warning(f"Could not parse rate constant {param_name}, using default 1.0")
                        k_values.append(1.0)
            
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
                possible_forward = [f'k{i+1}f', f'k{i+1}F', f'k{i+1}_f', f'k{i+1}', f'k{2*i+1}']
                possible_reverse = [f'k{i+1}r', f'k{i+1}R', f'k{i+1}_r', f'k{i+1}_rev', f'k{2*i+2}']
                
                kf_found = False
                kr_found = False
                
                # Look for forward rate constant
                for kf_name in possible_forward:
                    kf_row = params_df[params_df['Parameter'] == kf_name]
                    if not kf_row.empty:
                        try:
                            k_forward.append(float(kf_row['Values'].iloc[0]))
                            kf_found = True
                            break
                        except (ValueError, TypeError):
                            continue
                
                # Look for reverse rate constant
                for kr_name in possible_reverse:
                    kr_row = params_df[params_df['Parameter'] == kr_name]
                    if not kr_row.empty:
                        try:
                            k_reverse.append(float(kr_row['Values'].iloc[0]))
                            kr_found = True
                            break
                        except (ValueError, TypeError):
                            continue
                
                # Use defaults if we couldn't find them
                if not kf_found:
                    st.warning(f"Forward rate constant for reaction {i+1} not found, using default 1.0")
                    k_forward.append(1.0)
                    
                if not kr_found:
                    st.warning(f"Reverse rate constant for reaction {i+1} not found, using default 1.0")
                    k_reverse.append(1.0)
        
        # Make sure we have the right number of rate constants
        while len(k_forward) < len(sigma_values):
            k_forward.append(1.0)
        while len(k_reverse) < len(sigma_values):
            k_reverse.append(1.0)
            
        # Calculate individual equilibrium constants with error checking
        K_eq = []
        for i in range(len(sigma_values)):
            if k_reverse[i] != 0 and not np.isnan(k_reverse[i]) and not np.isinf(k_reverse[i]):
                K_i = k_forward[i] / k_reverse[i]
                if np.isnan(K_i) or np.isinf(K_i):
                    st.warning(f"Invalid equilibrium constant for reaction {i+1}, using 1.0")
                    K_i = 1.0
                K_eq.append(K_i)
            else:
                st.warning(f"Zero or invalid reverse rate constant for reaction {i+1}, using K_eq = 1.0")
                K_eq.append(1.0)
        
        # Calculate overall equilibrium constant: K_overall = ∏(K_i^σ_i)
        K_overall = 1.0
        for i, (K_i, sigma_i) in enumerate(zip(K_eq, sigma_values)):
            if K_i > 0 and not np.isnan(K_i) and not np.isinf(K_i):  # Make sure we can take the power
                try:
                    contribution = K_i ** sigma_i
                    if np.isnan(contribution) or np.isinf(contribution):
                        st.warning(f"Invalid contribution from reaction {i+1}, skipping")
                        continue
                    K_overall *= contribution
                except (OverflowError, ZeroDivisionError):
                    st.warning(f"Mathematical error in reaction {i+1}, using contribution = 1.0")
                    continue
            else:
                st.warning(f"Invalid equilibrium constant for reaction {i+1} in overall calculation")
        
        # Final validation of K_overall
        if np.isnan(K_overall) or np.isinf(K_overall):
            st.error("Overall equilibrium constant is invalid, using 1.0")
            K_overall = 1.0
        
        return K_eq, K_overall, T, R, k_forward, k_reverse
        
    except Exception as e:
        st.error(f"Error in equilibrium constant calculation: {str(e)}")
        # Return safe defaults if calculation fails
        n_reactions = len(sigma_values)
        return [1.0] * n_reactions, 1.0, 320.0, 8.314, [1.0] * n_reactions, [1.0] * n_reactions

def create_equilibrium_plots(K_eq_list, sigma_values, reaction_names, K_overall):
    """
    Creates comprehensive plots for equilibrium analysis
    """
    # Individual equilibrium constants bar chart
    fig_individual = go.Figure()
    
    # Color code based on equilibrium position
    colors = []
    for K_i in K_eq_list:
        if K_i > 10:
            colors.append('darkgreen')  # Strongly forward
        elif K_i > 1:
            colors.append('lightgreen')  # Forward
        elif K_i == 1:
            colors.append('yellow')  # Balanced
        elif K_i > 0.1:
            colors.append('orange')  # Reverse
        else:
            colors.append('red')  # Strongly reverse
    
    fig_individual.add_trace(go.Bar(
        x=reaction_names,
        y=K_eq_list,
        marker_color=colors,
        hovertemplate="<b>%{x}</b><br>K_eq: %{y:.2e}<br>σ: %{customdata:.3f}<extra></extra>",
        customdata=sigma_values,
        name="Individual K_eq"
    ))
    
    fig_individual.update_layout(
        title='Individual Equilibrium Constants by Reaction',
        xaxis_title='Reaction',
        yaxis_title='K_eq (log scale)',
        yaxis_type="log",
        template='plotly_white',
        height=500,
        showlegend=False
    )
    
    # Add horizontal line at K=1 (equilibrium)
    fig_individual.add_hline(y=1, line_dash="dash", line_color="black", 
                           annotation_text="K=1 (Equilibrium)")
    
    # Contribution to overall equilibrium (K_i^σ_i)
    contributions = [K_i**sigma_i for K_i, sigma_i in zip(K_eq_list, sigma_values)]
    
    fig_contributions = go.Figure()
    
    fig_contributions.add_trace(go.Bar(
        x=reaction_names,
        y=contributions,
        marker_color='steelblue',
        hovertemplate="<b>%{x}</b><br>K_eq^σ: %{y:.2e}<br>K_eq: %{customdata[0]:.2e}<br>σ: %{customdata[1]:.3f}<extra></extra>",
        customdata=list(zip(K_eq_list, sigma_values)),
        name="Contributions to K_overall"
    ))
    
    fig_contributions.update_layout(
        title='Contributions to Overall Equilibrium Constant (K_eq^σ)',
        xaxis_title='Reaction',
        yaxis_title='K_eq^σ (log scale)',
        yaxis_type="log",
        template='plotly_white',
        height=500,
        showlegend=False
    )
    
    return fig_individual, fig_contributions

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
        tab1, tab2, tab3, tab4, tab5 = st.tabs(["📋 Reaction Mechanism", "📊 Stoichiometric Matrix", "⚛️ Atomic Matrix", "⚙️ Parameters", "🌡️ Thermodynamics Analysis"])
        
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
            st.subheader("🌡️ Thermodynamics Analysis")
            
            # Explanation of what this analysis does
            st.markdown("""
            <div style="background-color: #f0f2f6; padding: 1rem; border-radius: 0.5rem; margin: 1rem 0; border-left: 4px solid #1f77b4;">
                <h4>🔬 Comprehensive Equilibrium Constant Analysis</h4>
                <p>This analysis calculates thermodynamic equilibrium constants from your kinetic parameters:</p>
                <ul>
                    <li><strong>Individual K_eq:</strong> For each reaction: K_i = k_forward/k_reverse</li>
                    <li><strong>Overall K_eq:</strong> Combined equilibrium: K_overall = ∏K_i^σ_i</li>
                    <li><strong>Thermodynamic insight:</strong> Understand which direction is thermodynamically favored</li>
                    <li><strong>Visual analysis:</strong> Interactive plots and comprehensive interpretation</li>
                </ul>
            </div>
            """, unsafe_allow_html=True)
            
            if st.session_state.param_data:
                # Convert our parameter data to a DataFrame for processing
                params_df = pd.DataFrame(st.session_state.param_data)
                
                # Parameter inspection section
                with st.expander("🔍 Parameter Inspection"):
                    col1, col2 = st.columns(2)
                    
                    with col1:
                        st.write("**System Parameters:**")
                        # Temperature and gas constant
                        temp_row = params_df[params_df['Parameter'].str.lower() == 't']
                        if not temp_row.empty:
                            temp_val = float(temp_row['Values'].iloc[0])
                            st.metric("Temperature", f"{temp_val:.1f} K")
                        else:
                            st.warning("Temperature (T) not found")
                        
                        r_row = params_df[params_df['Parameter'].str.lower() == 'r']
                        if not r_row.empty:
                            r_val = float(r_row['Values'].iloc[0])
                            st.metric("Gas Constant", f"{r_val:.3f} J/(mol·K)")
                        else:
                            st.warning("Gas constant (R) not found")
                    
                    with col2:
                        st.write("**Rate Constants Analysis:**")
                        # Analyze rate constants
                        k_params = params_df[params_df['Parameter'].str.contains(r'^k\d+', regex=True, na=False)]
                        if not k_params.empty:
                            st.metric("Rate Constants Found", len(k_params))
                            
                            # Check for forward/reverse pattern
                            forward_k = k_params[k_params['Parameter'].str.endswith('f')]
                            reverse_k = k_params[k_params['Parameter'].str.endswith('r')]
                            
                            if len(forward_k) > 0 and len(reverse_k) > 0:
                                st.success("✅ Forward/Reverse pattern detected")
                            else:
                                st.info("ℹ️ Using alternating pattern detection")
                        else:
                            st.error("❌ No rate constants found")
                    
                    # Show rate constants table
                    if not k_params.empty:
                        st.write("**All Rate Constants:**")
                        st.dataframe(k_params[['Parameter', 'Values', 'Units']], use_container_width=True)
                
                # Rate constant pattern selection
                st.markdown("### 🔧 Rate Constant Configuration")
                
                col1, col2 = st.columns(2)
                
                with col1:
                    k_pattern_gen = st.radio(
                        "Rate constant organization:",
                        ["auto", "alternating"],
                        format_func=lambda x: {
                            "auto": "Auto-detect (kXf/kXr pattern)",
                            "alternating": "Alternating (k1=forward, k2=reverse, k3=forward, k4=reverse...)"
                        }[x],
                        help="Select how your rate constants are named/organized",
                        key="gen_k_pattern"
                    )
                
                with col2:
                    # Show detected pattern
                    k_params = params_df[params_df['Parameter'].str.contains(r'k', regex=True, na=False)]
                    if not k_params.empty:
                        k_names = k_params['Parameter'].tolist()
                        has_f_suffix = any('f' in k.lower() for k in k_names)
                        has_r_suffix = any('r' in k.lower() for k in k_names)
                        
                        if has_f_suffix and has_r_suffix:
                            st.success("✅ Detected: Forward/Reverse suffixes")
                        else:
                            st.info("ℹ️ Detected: Sequential numbering")
                
                # Parameter editing with validation
                st.markdown("### ✏️ Edit Parameters")
                
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
                
                # Stoichiometric numbers section
                st.markdown("---")
                st.subheader("📏 Stoichiometric Numbers (σ)")
                
                # Explanation
                st.markdown("""
                **What are stoichiometric numbers?**
                These coefficients determine how each elementary reaction contributes to the overall process:
                - **σ = 1**: Reaction occurs once per overall cycle
                - **σ = 2**: Reaction occurs twice per overall cycle  
                - **σ = 0**: Reaction is in side equilibrium
                - **σ = -1**: Reaction runs in reverse direction
                """)
                
                # Get the number of reactions for sigma input
                num_reactions = len(st.session_state.reactions)
                
                # Show current reactions for reference
                with st.expander("📋 Current Reaction Mechanism"):
                    for i, (_, reactants, products, original) in enumerate(st.session_state.reactions):
                        st.write(f"**r{i+1}:** {original}")
                
                # Initialize sigma values
                if f'sigma_values_gen' not in st.session_state:
                    st.session_state.sigma_values_gen = [1.0] * num_reactions
                elif len(st.session_state.sigma_values_gen) != num_reactions:
                    st.session_state.sigma_values_gen = [1.0] * num_reactions
                
                # Create input fields for sigma values
                st.write("**Enter stoichiometric number for each reaction:**")
                sigma_inputs_gen = []
                
                # Create a more organized layout
                for i in range(num_reactions):
                    col1, col2, col3 = st.columns([2, 1, 3])
                    
                    with col1:
                        original_reaction = st.session_state.reactions[i][3] if i < len(st.session_state.reactions) else f"r{i+1}"
                        st.write(f"**r{i+1}:**")
                    
                    with col2:
                        sigma_val = st.number_input(
                            f"σ{i+1}",
                            value=st.session_state.sigma_values_gen[i],
                            step=0.1,
                            format="%.3f",
                            key=f"gen_sigma_{i}",
                            help=f"Stoichiometric number for reaction r{i+1}"
                        )
                        sigma_inputs_gen.append(sigma_val)
                    
                    with col3:
                        # Show the actual reaction
                        st.code(original_reaction, language=None)
                
                st.session_state.sigma_values_gen = sigma_inputs_gen
                
                # Validation section
                st.markdown("### ✅ Pre-calculation Validation")
                
                col1, col2, col3 = st.columns(3)
                
                with col1:
                    # Rate constant validation
                    k_params = edited_params_gen[edited_params_gen['Parameter'].str.contains(r'k', regex=True, na=False)]
                    expected_k_count = num_reactions * 2  # Forward and reverse for each reaction
                    
                    if len(k_params) >= expected_k_count:
                        st.success(f"✅ Rate constants: {len(k_params)}/{expected_k_count}")
                    else:
                        st.warning(f"⚠️ Rate constants: {len(k_params)}/{expected_k_count}")
                
                with col2:
                    # Temperature validation
                    temp_params = edited_params_gen[edited_params_gen['Parameter'].str.lower() == 't']
                    if len(temp_params) > 0:
                        temp_val = float(temp_params['Values'].iloc[0])
                        if 200 <= temp_val <= 2000:
                            st.success(f"✅ Temperature: {temp_val:.1f} K")
                        else:
                            st.warning(f"⚠️ Temperature: {temp_val:.1f} K (unusual)")
                    else:
                        st.error("❌ Temperature missing")
                
                with col3:
                    # Sigma validation
                    sigma_sum = sum(st.session_state.sigma_values_gen)
                    st.info(f"ℹ️ Σσ = {sigma_sum:.3f}")
                    if abs(sigma_sum) < 0.1:
                        st.warning("Very small Σσ - check if correct")
                
                # Main calculation button
                if st.button("🧮 Calculate Equilibrium Analysis", type="primary", key="gen_calc_eq"):
                    
                    with st.spinner("Performing equilibrium analysis..."):
                        try:
                            # Calculate equilibrium constants using function
                            K_eq_list_gen, K_overall_gen, T_gen, R_gen, k_forward_gen, k_reverse_gen = calculate_equilibrium_constants(
                                edited_params_gen, st.session_state.sigma_values_gen, k_pattern_gen
                            )
                            
                            # Store comprehensive results
                            st.session_state.equilibrium_results = {
                                'K_eq_list': K_eq_list_gen,
                                'K_overall': K_overall_gen,
                                'sigma_values': st.session_state.sigma_values_gen,
                                'k_forward': k_forward_gen,
                                'k_reverse': k_reverse_gen,
                                'temperature': T_gen,
                                'gas_constant': R_gen,
                                'reaction_names': [f"r{i+1}" for i in range(len(K_eq_list_gen))]
                            }
                            
                            st.success("✅ Equilibrium analysis completed!")
                            
                            # Results display
                            st.markdown("---")
                            st.markdown("## 🎯 Equilibrium Analysis Results")
                            
                            # Summary metrics at the top
                            col1, col2, col3, col4 = st.columns(4)
                            
                            with col1:
                                st.metric("Temperature", f"{T_gen:.1f} K")
                            
                            with col2:
                                if K_overall_gen > 0 and not np.isinf(K_overall_gen):
                                    ln_K = math.log(K_overall_gen)
                                    st.metric("ln(K_overall)", f"{ln_K:.4f}")
                                else:
                                    st.metric("ln(K_overall)", "Undefined")
                            
                            with col3:
                                if K_overall_gen > 0 and not np.isinf(K_overall_gen):
                                    # Calculate standard Gibbs free energy change
                                    delta_G = -R_gen * T_gen * math.log(K_overall_gen) / 1000  # Convert to kJ/mol
                                    st.metric("ΔG° (kJ/mol)", f"{delta_G:.2f}")
                                else:
                                    st.metric("ΔG° (kJ/mol)", "Undefined")
                            
                            with col4:
                                # Show overall equilibrium constant with scientific notation
                                if K_overall_gen != 0 and not np.isinf(K_overall_gen):
                                    exponent = math.floor(math.log10(abs(K_overall_gen)))
                                    mantissa = K_overall_gen / (10 ** exponent)
                                    st.metric("K_overall", f"{mantissa:.2f}×10^{exponent}")
                                else:
                                    st.metric("K_overall", f"{K_overall_gen}")
                            
                            # Thermodynamic interpretation
                            st.subheader("🔬 Thermodynamic Interpretation")
                            
                            col1, col2 = st.columns(2)
                            
                            with col1:
                                st.markdown("**Overall Process:**")
                                if K_overall_gen > 1000:
                                    st.success("⚡ **Strongly favors products** - Excellent thermodynamics")
                                elif K_overall_gen > 1:
                                    st.info("➡️ **Favors products** - Good thermodynamics")
                                elif K_overall_gen == 1:
                                    st.warning("⚖️ **Equilibrium balanced** - Marginal thermodynamics")
                                elif K_overall_gen > 0.001:
                                    st.info("⬅️ **Favors reactants** - May need optimization")
                                else:
                                    st.error("⚡ **Strongly favors reactants** - Poor thermodynamics")
                            
                            with col2:
                                st.markdown("**Process Feasibility:**")
                                if K_overall_gen > 100:
                                    st.success("✅ **Highly feasible** - Should proceed well")
                                elif K_overall_gen > 1:
                                    st.success("✅ **Feasible** - Favorable conditions")
                                elif K_overall_gen > 0.01:
                                    st.warning("⚠️ **Marginally feasible** - Consider optimization")
                                else:
                                    st.error("❌ **Not feasible** - Requires extreme conditions")
                            
                            # Individual reaction analysis table
                            st.subheader("📊 Individual Reaction Analysis")
                            
                            eq_data = []
                            for i, (K_i, sigma_i) in enumerate(zip(K_eq_list_gen, st.session_state.sigma_values_gen)):
                                # Get the original reaction for reference
                                original_reaction = st.session_state.reactions[i][3] if i < len(st.session_state.reactions) else f"r{i+1}"
                                
                                # Determine thermodynamic favorability
                                if K_i > 100:
                                    favorability = "Strongly Forward"
                                    color = "🟢"
                                elif K_i > 1:
                                    favorability = "Forward"
                                    color = "🔵"
                                elif K_i == 1:
                                    favorability = "Balanced"
                                    color = "🟡"
                                elif K_i > 0.01:
                                    favorability = "Reverse"
                                    color = "🟠"
                                else:
                                    favorability = "Strongly Reverse"
                                    color = "🔴"
                                
                                eq_data.append({
                                    'Reaction': f"r{i+1}",
                                    'Original_Equation': original_reaction,
                                    'k_forward': f'{k_forward_gen[i]:.2e}',
                                    'k_reverse': f'{k_reverse_gen[i]:.2e}',
                                    'K_eq': f'{K_i:.2e}',
                                    'σ': f'{sigma_i:.3f}',
                                    'K_eq^σ': f'{K_i**sigma_i:.2e}',
                                    'Favorability': f'{color} {favorability}'
                                })
                            
                            eq_df = pd.DataFrame(eq_data)
                            st.dataframe(eq_df, use_container_width=True)
                            
                            # Visual analysis
                            st.subheader("📈 Visual Analysis")
                            
                            # Get reaction names for plotting
                            reaction_names = [f"r{i+1}: {st.session_state.reactions[i][3][:20]}..." if len(st.session_state.reactions[i][3]) > 20 
                                            else f"r{i+1}: {st.session_state.reactions[i][3]}" 
                                            for i in range(len(K_eq_list_gen))]
                            
                            # Create plots
                            fig_individual, fig_contributions = create_equilibrium_plots(
                                K_eq_list_gen, st.session_state.sigma_values_gen, 
                                reaction_names, K_overall_gen
                            )
                            
                            # Display plots
                            st.plotly_chart(fig_individual, use_container_width=True)
                            st.plotly_chart(fig_contributions, use_container_width=True)
                            
                            # Mechanism analysis
                            st.subheader("🔬 Mechanism Analysis")
                            
                            # Find rate-limiting and driving steps
                            min_K = min(K_eq_list_gen)
                            min_K_idx = K_eq_list_gen.index(min_K)
                            min_K_reaction = f"r{min_K_idx+1}"
                            
                            max_K = max(K_eq_list_gen)
                            max_K_idx = K_eq_list_gen.index(max_K)
                            max_K_reaction = f"r{max_K_idx+1}"
                            
                            col1, col2 = st.columns(2)
                            
                            with col1:
                                st.markdown("**🚧 Potential Thermodynamic Bottlenecks:**")
                                if min_K < 0.1:
                                    st.write(f"• **{min_K_reaction}**: K = {min_K:.2e} (thermodynamically unfavorable)")
                                    st.write(f"  Original: {st.session_state.reactions[min_K_idx][3]}")
                                    st.write("  💡 **Recommendation:** Optimize conditions for this step")
                                else:
                                    st.write("• ✅ No major thermodynamic bottlenecks detected")
                                
                                # Check for highly unfavorable steps
                                unfavorable_steps = [i for i, K in enumerate(K_eq_list_gen) if K < 0.01]
                                if unfavorable_steps:
                                    st.write("**Additional unfavorable steps:**")
                                    for idx in unfavorable_steps:
                                        st.write(f"• r{idx+1}: K = {K_eq_list_gen[idx]:.2e}")
                            
                            with col2:
                                st.markdown("**⚡ Thermodynamic Driving Forces:**")
                                if max_K > 10:
                                    st.write(f"• **{max_K_reaction}**: K = {max_K:.2e} (strongly favored)")
                                    st.write(f"  Original: {st.session_state.reactions[max_K_idx][3]}")
                                    st.write("  ✅ **This step provides strong driving force**")
                                else:
                                    st.write("• All steps have moderate equilibrium constants")
                                
                                # Check for other highly favorable steps
                                favorable_steps = [i for i, K in enumerate(K_eq_list_gen) if K > 100]
                                if len(favorable_steps) > 1:
                                    st.write("**Additional driving forces:**")
                                    for idx in favorable_steps:
                                        if idx != max_K_idx:
                                            st.write(f"• r{idx+1}: K = {K_eq_list_gen[idx]:.2e}")
                            
                            # Practical recommendations
                            st.subheader("💡 Practical Recommendations")
                            
                            recommendations = []
                            
                            # Overall thermodynamics recommendations
                            if K_overall_gen > 1000:
                                recommendations.append("✅ **Excellent thermodynamics** - Reaction should proceed well under standard conditions")
                                recommendations.append("🎯 **Focus on kinetics** - Thermodynamics are favorable, optimize rate constants and mass transfer")
                            elif K_overall_gen > 1:
                                recommendations.append("✅ **Good thermodynamics** - Reaction is thermodynamically favorable")
                                recommendations.append("⚖️ **Balanced optimization** - Consider both thermodynamic and kinetic improvements")
                            elif K_overall_gen > 0.001:
                                recommendations.append("⚠️ **Marginal thermodynamics** - Consider temperature, pressure, or concentration optimization")
                                recommendations.append("🌡️ **Temperature sensitivity** - Higher temperatures might improve equilibrium position")
                            else:
                                recommendations.append("❌ **Poor thermodynamics** - Reaction may require extreme conditions or process coupling")
                                recommendations.append("🔄 **Process redesign** - Consider alternative pathways or reaction coupling")
                            
                            # Specific step recommendations
                            unfavorable_steps = [i for i, K in enumerate(K_eq_list_gen) if K < 0.01]
                            if unfavorable_steps:
                                step_names = [f"r{i+1}" for i in unfavorable_steps]
                                recommendations.append(f"🔧 **Optimize specific steps:** {', '.join(step_names)} have very low K_eq values")
                                recommendations.append("💡 **Consider:** Alternative catalysts, modified conditions, or step elimination")
                            
                            # Temperature optimization suggestions
                            if any(K > 1000 for K in K_eq_list_gen) and any(K < 0.001 for K in K_eq_list_gen):
                                recommendations.append("🌡️ **Temperature optimization** may help balance competing thermodynamic requirements")
                                recommendations.append("📊 **Consider:** Van 't Hoff analysis to determine optimal temperature")
                            
                            # Process integration recommendations
                            if 0.1 < K_overall_gen < 10:
                                recommendations.append("🔗 **Process integration** - Consider coupling with thermodynamically favorable reactions")
                                recommendations.append("♻️ **Separation strategy** - Product removal might shift equilibrium favorably")
                            
                            # Display recommendations
                            for i, rec in enumerate(recommendations, 1):
                                st.write(f"{i}. {rec}")
                            
                            # Download section
                            st.subheader("💾 Download Results")
                            
                            # Create comprehensive results for download
                            download_data = []
                            for i in range(len(K_eq_list_gen)):
                                original_reaction = st.session_state.reactions[i][3] if i < len(st.session_state.reactions) else f"r{i+1}"
                                download_data.append({
                                    'Reaction_ID': f"r{i+1}",
                                    'Original_Equation': original_reaction,
                                    'k_forward': k_forward_gen[i],
                                    'k_reverse': k_reverse_gen[i],
                                    'K_equilibrium': K_eq_list_gen[i],
                                    'Sigma': st.session_state.sigma_values_gen[i],
                                    'K_eq_to_sigma': K_eq_list_gen[i]**st.session_state.sigma_values_gen[i],
                                    'ln_K_eq': math.log(K_eq_list_gen[i]) if K_eq_list_gen[i] > 0 else 'Undefined',
                                    'Delta_G_individual_kJ_mol': -R_gen * T_gen * math.log(K_eq_list_gen[i]) / 1000 if K_eq_list_gen[i] > 0 else 'Undefined',
                                    'Thermodynamic_Favorability': 'Forward' if K_eq_list_gen[i] > 1 else 'Reverse' if K_eq_list_gen[i] < 1 else 'Balanced'
                                })
                            
                            # Add overall summary
                            overall_summary = {
                                'Reaction_ID': 'OVERALL',
                                'Original_Equation': 'Combined Process',
                                'k_forward': '',
                                'k_reverse': '',
                                'K_equilibrium': K_overall_gen,
                                'Sigma': '',
                                'K_eq_to_sigma': '',
                                'ln_K_eq': math.log(K_overall_gen) if K_overall_gen > 0 else 'Undefined',
                                'Delta_G_individual_kJ_mol': -R_gen * T_gen * math.log(K_overall_gen) / 1000 if K_overall_gen > 0 else 'Undefined',
                                'Thermodynamic_Favorability': 'Forward' if K_overall_gen > 1 else 'Reverse' if K_overall_gen < 1 else 'Balanced'
                            }
                            download_data.append(overall_summary)
                            
                            # Add system parameters
                            system_params = {
                                'Reaction_ID': 'SYSTEM_PARAMS',
                                'Original_Equation': f'T={T_gen}K, R={R_gen}J/(mol·K)',
                                'k_forward': f'Sum_sigma={sum(st.session_state.sigma_values_gen):.3f}',
                                'k_reverse': f'Num_reactions={len(K_eq_list_gen)}',
                                'K_equilibrium': '',
                                'Sigma': '',
                                'K_eq_to_sigma': '',
                                'ln_K_eq': '',
                                'Delta_G_individual_kJ_mol': '',
                                'Thermodynamic_Favorability': ''
                            }
                            download_data.append(system_params)
                            
                            results_df = pd.DataFrame(download_data)
                            csv = results_df.to_csv(index=False)
                            
                            st.download_button(
                                label="📄 Download Equilibrium Analysis CSV",
                                data=csv,
                                file_name="thermodynamic_equilibrium_analysis.csv",
                                mime="text/csv",
                                help="Download comprehensive equilibrium analysis results with thermodynamic properties"
                            )
                            
                        except Exception as e:
                            st.error(f"❌ Error in equilibrium analysis: {str(e)}")
                            with st.expander("🔧 Error Details"):
                                st.write(str(e))
                                st.write("**Troubleshooting:**")
                                st.write("• Verify all rate constants are positive numbers")
                                st.write("• Check that forward and reverse rate constants are properly paired")
                                st.write("• Ensure temperature is reasonable (200-2000 K)")
                                st.write("• Validate stoichiometric numbers are appropriate")
                
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

    else:
        # Welcome screen when no reactions have been added yet
        st.markdown("""
        ## 🚀 Getting Started
        
        Enter chemical reactions **one at a time** in the sidebar to automatically generate:
        - **Stoichiometric Matrix**: Shows how each reaction affects each species
        - **Atomic Matrix**: Shows atomic composition of each species  
        - **Parameter Table**: Editable parameters for reaction kinetics
        - **Thermodynamics Analysis**: Comprehensive equilibrium constant calculations with visual analysis
        
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
        - **Thermodynamics**: Calculate K_overall = ∏K_i^σ_i with comprehensive analysis
        - **Visual Analysis**: Interactive plots and detailed interpretations
        - **Smart Recommendations**: Get practical suggestions for your mechanism
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
        - **Thermodynamics**: Comprehensive equilibrium constant analysis with K_overall = ∏K_i^σ_i
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
                
                # Explanation of what this analysis does
                st.markdown("""
                <div style="background-color: #f0f2f6; padding: 1rem; border-radius: 0.5rem; margin: 1rem 0; border-left: 4px solid #1f77b4;">
                    <h4>🔬 Comprehensive Equilibrium Constant Analysis</h4>
                    <p>This analysis calculates thermodynamic equilibrium constants from your kinetic parameters:</p>
                    <ul>
                        <li><strong>Individual K_eq:</strong> For each reaction: K_i = k_forward/k_reverse</li>
                        <li><strong>Overall K_eq:</strong> Combined equilibrium: K_overall = ∏K_i^σ_i</li>
                        <li><strong>Thermodynamic insight:</strong> Understand which direction is thermodynamically favored</li>
                        <li><strong>Visual analysis:</strong> Interactive plots and comprehensive interpretation</li>
                    </ul>
                </div>
                """, unsafe_allow_html=True)
                
                if 'params' in loaded_files:
                    # Display and allow editing of parameters
                    st.subheader("Parameters")
                    
                    # Parameter inspection section
                    with st.expander("🔍 Parameter Inspection"):
                        param_preview = loaded_files['params'].copy()
                        
                        col1, col2 = st.columns(2)
                        
                        with col1:
                            st.write("**System Parameters:**")
                            # Temperature and gas constant
                            temp_row = param_preview[param_preview['Parameter'].str.lower() == 't']
                            if not temp_row.empty:
                                temp_val = float(temp_row['Values'].iloc[0])
                                st.metric("Temperature", f"{temp_val:.1f} K")
                            else:
                                st.warning("Temperature (T) not found")
                            
                            r_row = param_preview[param_preview['Parameter'].str.lower() == 'r']
                            if not r_row.empty:
                                r_val = float(r_row['Values'].iloc[0])
                                st.metric("Gas Constant", f"{r_val:.3f} J/(mol·K)")
                            else:
                                st.warning("Gas constant (R) not found")
                        
                        with col2:
                            st.write("**Rate Constants Analysis:**")
                            # Analyze rate constants
                            k_params = param_preview[param_preview['Parameter'].str.contains(r'k', regex=True, na=False)]
                            if not k_params.empty:
                                st.metric("Rate Constants Found", len(k_params))
                                
                                # Check for forward/reverse pattern
                                forward_k = k_params[k_params['Parameter'].str.endswith('f')]
                                reverse_k = k_params[k_params['Parameter'].str.endswith('r')]
                                
                                if len(forward_k) > 0 and len(reverse_k) > 0:
                                    st.success("✅ Forward/Reverse pattern detected")
                                    st.write(f"Forward: {len(forward_k)}, Reverse: {len(reverse_k)}")
                                else:
                                    st.info("ℹ️ Using alternating pattern detection")
                            else:
                                st.error("❌ No rate constants found")
                        
                        # Show all parameters table
                        st.write("**All Parameters:**")
                        st.dataframe(param_preview, use_container_width=True)
                    
                    # Rate constant pattern selection
                    st.markdown("### 🔧 Rate Constant Configuration")
                    
                    k_pattern = st.radio(
                        "Rate constant organization:",
                        ["auto", "alternating"],
                        format_func=lambda x: {
                            "auto": "Auto-detect (kXf/kXr or similar)",
                            "alternating": "Alternating (k1=forward, k2=reverse, k3=forward, k4=reverse...)"
                        }[x],
                        help="Select how your rate constants are named/organized"
                    )
                    
                    # Editable parameters table
                    st.markdown("### ✏️ Edit Parameters")
                    
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
                        st.subheader("📏 Stoichiometric Numbers (σ)")
                        st.markdown("Enter the stoichiometric number for each reaction step:")
                        
                        # Explanation
                        with st.expander("📖 Understanding Stoichiometric Numbers"):
                            st.markdown("""
                            **Stoichiometric numbers (σ) determine how each elementary reaction contributes to the overall process:**
                            
                            - **σ = 1**: Reaction occurs once per overall cycle
                            - **σ = 2**: Reaction occurs twice per overall cycle  
                            - **σ = 0**: Reaction is in side equilibrium (doesn't contribute)
                            - **σ = -1**: Reaction runs in reverse direction
                            
                            **For a simple linear mechanism:** σ = [1, 1, 1, ...] (all steps occur once)
                            **For branched mechanisms:** Some σ values may be different
                            """)
                        
                        # Initialize sigma values if needed
                        if len(st.session_state.sigma_values) != num_reactions:
                            st.session_state.sigma_values = [1.0] * num_reactions
                        
                        # Create input fields for sigma values
                        sigma_inputs = []
                        
                        # Show reaction info if available
                        if 'stoich' in loaded_files and len(loaded_files['stoich']) > 0:
                            with st.expander("📋 Reaction Overview"):
                                for i in range(num_reactions):
                                    if i < len(loaded_files['stoich']):
                                        reaction_name = loaded_files['stoich'].iloc[i, 0]
                                        st.write(f"**r{i+1}:** {reaction_name}")
                        
                        # Organized sigma input
                        cols = st.columns(min(4, num_reactions))
                        for i in range(num_reactions):
                            col_idx = i % 4
                            with cols[col_idx]:
                                reaction_name = loaded_files['stoich'].iloc[i, 0] if i < len(loaded_files['stoich']) else f"r{i+1}"
                                sigma_val = st.number_input(
                                    f"σ{i+1} ({reaction_name})",
                                    value=st.session_state.sigma_values[i],
                                    step=0.1,
                                    format="%.3f",
                                    key=f"analysis_sigma_{i}",
                                    help=f"Stoichiometric number for {reaction_name}"
                                )
                                sigma_inputs.append(sigma_val)
                        
                        st.session_state.sigma_values = sigma_inputs
                        
                        # Validation section
                        st.markdown("### ✅ Pre-calculation Validation")
                        
                        col1, col2, col3 = st.columns(3)
                        
                        with col1:
                            # Rate constant validation
                            k_params = edited_params[edited_params['Parameter'].str.contains(r'k', regex=True, na=False)]
                            expected_k_count = num_reactions * 2  # Forward and reverse for each reaction
                            
                            if len(k_params) >= expected_k_count:
                                st.success(f"✅ Rate constants: {len(k_params)}/{expected_k_count}")
                            else:
                                st.warning(f"⚠️ Rate constants: {len(k_params)}/{expected_k_count}")
                        
                        with col2:
                            # Temperature validation
                            temp_params = edited_params[edited_params['Parameter'].str.lower() == 't']
                            if len(temp_params) > 0:
                                temp_val = float(temp_params['Values'].iloc[0])
                                if 200 <= temp_val <= 2000:
                                    st.success(f"✅ Temperature: {temp_val:.1f} K")
                                else:
                                    st.warning(f"⚠️ Temperature: {temp_val:.1f} K (unusual)")
                            else:
                                st.error("❌ Temperature missing")
                        
                        with col3:
                            # Sigma validation
                            sigma_sum = sum(st.session_state.sigma_values)
                            st.info(f"ℹ️ Σσ = {sigma_sum:.3f}")
                            if abs(sigma_sum) < 0.1:
                                st.warning("Very small Σσ - check if correct")
                        
                        # Calculation button
                        if st.button("🧮 Calculate Equilibrium Analysis", type="primary"):
                            
                            with st.spinner("Performing equilibrium analysis..."):
                                try:
                                    # Use calculation function
                                    K_eq_list, K_overall, T, R, k_forward, k_reverse = calculate_equilibrium_constants(
                                        edited_params, st.session_state.sigma_values, k_pattern
                                    )
                                    
                                    # Store comprehensive results
                                    st.session_state.equilibrium_results = {
                                        'K_eq_list': K_eq_list,
                                        'K_overall': K_overall,
                                        'sigma_values': st.session_state.sigma_values,
                                        'k_forward': k_forward,
                                        'k_reverse': k_reverse,
                                        'temperature': T,
                                        'gas_constant': R,
                                        'reaction_names': [f"r{i+1}" for i in range(len(K_eq_list))]
                                    }
                                    
                                    st.success("✅ Equilibrium analysis completed!")
                                    
                                    # Results display - complete analysis
                                    st.markdown("---")
                                    st.markdown("## 🎯 Equilibrium Analysis Results")
                                    
                                    # Summary metrics at the top
                                    col1, col2, col3, col4 = st.columns(4)
                                    
                                    with col1:
                                        st.metric("Temperature", f"{T:.1f} K")
                                    
                                    with col2:
                                        if K_overall > 0 and not np.isinf(K_overall):
                                            ln_K = math.log(K_overall)
                                            st.metric("ln(K_overall)", f"{ln_K:.4f}")
                                        else:
                                            st.metric("ln(K_overall)", "Undefined")
                                    
                                    with col3:
                                        if K_overall > 0 and not np.isinf(K_overall):
                                            # Calculate standard Gibbs free energy change
                                            delta_G = -R * T * math.log(K_overall) / 1000  # Convert to kJ/mol
                                            st.metric("ΔG° (kJ/mol)", f"{delta_G:.2f}")
                                        else:
                                            st.metric("ΔG° (kJ/mol)", "Undefined")
                                    
                                    with col4:
                                        # Show overall equilibrium constant with scientific notation
                                        if K_overall != 0 and not np.isinf(K_overall):
                                            exponent = math.floor(math.log10(abs(K_overall)))
                                            mantissa = K_overall / (10 ** exponent)
                                            st.metric("K_overall", f"{mantissa:.2f}×10^{exponent}")
                                        else:
                                            st.metric("K_overall", f"{K_overall}")
                                    
                                    # Thermodynamic interpretation
                                    st.subheader("🔬 Thermodynamic Interpretation")
                                    
                                    col1, col2 = st.columns(2)
                                    
                                    with col1:
                                        st.markdown("**Overall Process:**")
                                        if K_overall > 1000:
                                            st.success("⚡ **Strongly favors products** - Excellent thermodynamics")
                                        elif K_overall > 1:
                                            st.info("➡️ **Favors products** - Good thermodynamics")
                                        elif K_overall == 1:
                                            st.warning("⚖️ **Equilibrium balanced** - Marginal thermodynamics")
                                        elif K_overall > 0.001:
                                            st.info("⬅️ **Favors reactants** - May need optimization")
                                        else:
                                            st.error("⚡ **Strongly favors reactants** - Poor thermodynamics")
                                    
                                    with col2:
                                        st.markdown("**Process Feasibility:**")
                                        if K_overall > 100:
                                            st.success("✅ **Highly feasible** - Should proceed well")
                                        elif K_overall > 1:
                                            st.success("✅ **Feasible** - Favorable conditions")
                                        elif K_overall > 0.01:
                                            st.warning("⚠️ **Marginally feasible** - Consider optimization")
                                        else:
                                            st.error("❌ **Not feasible** - Requires extreme conditions")
                                    
                                    # Individual reaction analysis table
                                    st.subheader("📊 Individual Reaction Analysis")
                                    
                                    eq_data = []
                                    for i, (K_i, sigma_i) in enumerate(zip(K_eq_list, st.session_state.sigma_values)):
                                        # Get the reaction name from stoich matrix if available
                                        if 'stoich' in loaded_files and i < len(loaded_files['stoich']):
                                            reaction_name = loaded_files['stoich'].iloc[i, 0]
                                        else:
                                            reaction_name = f"r{i+1}"
                                        
                                        # Determine thermodynamic favorability
                                        if K_i > 100:
                                            favorability = "Strongly Forward"
                                            color = "🟢"
                                        elif K_i > 1:
                                            favorability = "Forward"
                                            color = "🔵"
                                        elif K_i == 1:
                                            favorability = "Balanced"
                                            color = "🟡"
                                        elif K_i > 0.01:
                                            favorability = "Reverse"
                                            color = "🟠"
                                        else:
                                            favorability = "Strongly Reverse"
                                            color = "🔴"
                                        
                                        eq_data.append({
                                            'Reaction': f"r{i+1}",
                                            'Reaction_Name': reaction_name,
                                            'k_forward': f'{k_forward[i]:.2e}',
                                            'k_reverse': f'{k_reverse[i]:.2e}',
                                            'K_eq': f'{K_i:.2e}',
                                            'σ': f'{sigma_i:.3f}',
                                            'K_eq^σ': f'{K_i**sigma_i:.2e}',
                                            'Favorability': f'{color} {favorability}'
                                        })
                                    
                                    eq_df = pd.DataFrame(eq_data)
                                    st.dataframe(eq_df, use_container_width=True)
                                    
                                    # Visual analysis
                                    st.subheader("📈 Visual Analysis")
                                    
                                    # Get reaction names for plotting
                                    reaction_names = []
                                    for i in range(len(K_eq_list)):
                                        if 'stoich' in loaded_files and i < len(loaded_files['stoich']):
                                            reaction_name = loaded_files['stoich'].iloc[i, 0]
                                            if len(reaction_name) > 20:
                                                reaction_names.append(f"r{i+1}: {reaction_name[:20]}...")
                                            else:
                                                reaction_names.append(f"r{i+1}: {reaction_name}")
                                        else:
                                            reaction_names.append(f"r{i+1}")
                                    
                                    # Create plots
                                    fig_individual, fig_contributions = create_equilibrium_plots(
                                        K_eq_list, st.session_state.sigma_values, 
                                        reaction_names, K_overall
                                    )
                                    
                                    # Display plots
                                    st.plotly_chart(fig_individual, use_container_width=True)
                                    st.plotly_chart(fig_contributions, use_container_width=True)
                                    
                                    # Mechanism analysis
                                    st.subheader("🔬 Mechanism Analysis")
                                    
                                    # Find rate-limiting and driving steps
                                    min_K = min(K_eq_list)
                                    min_K_idx = K_eq_list.index(min_K)
                                    min_K_reaction = f"r{min_K_idx+1}"
                                    
                                    max_K = max(K_eq_list)
                                    max_K_idx = K_eq_list.index(max_K)
                                    max_K_reaction = f"r{max_K_idx+1}"
                                    
                                    col1, col2 = st.columns(2)
                                    
                                    with col1:
                                        st.markdown("**🚧 Potential Thermodynamic Bottlenecks:**")
                                        if min_K < 0.1:
                                            reaction_name = loaded_files['stoich'].iloc[min_K_idx, 0] if 'stoich' in loaded_files and min_K_idx < len(loaded_files['stoich']) else f"r{min_K_idx+1}"
                                            st.write(f"• **{min_K_reaction}**: K = {min_K:.2e} (thermodynamically unfavorable)")
                                            st.write(f"  Reaction: {reaction_name}")
                                            st.write("  💡 **Recommendation:** Optimize conditions for this step")
                                        else:
                                            st.write("• ✅ No major thermodynamic bottlenecks detected")
                                        
                                        # Check for highly unfavorable steps
                                        unfavorable_steps = [i for i, K in enumerate(K_eq_list) if K < 0.01]
                                        if unfavorable_steps:
                                            st.write("**Additional unfavorable steps:**")
                                            for idx in unfavorable_steps:
                                                st.write(f"• r{idx+1}: K = {K_eq_list[idx]:.2e}")
                                    
                                    with col2:
                                        st.markdown("**⚡ Thermodynamic Driving Forces:**")
                                        if max_K > 10:
                                            reaction_name = loaded_files['stoich'].iloc[max_K_idx, 0] if 'stoich' in loaded_files and max_K_idx < len(loaded_files['stoich']) else f"r{max_K_idx+1}"
                                            st.write(f"• **{max_K_reaction}**: K = {max_K:.2e} (strongly favored)")
                                            st.write(f"  Reaction: {reaction_name}")
                                            st.write("  ✅ **This step provides strong driving force**")
                                        else:
                                            st.write("• All steps have moderate equilibrium constants")
                                        
                                        # Check for other highly favorable steps
                                        favorable_steps = [i for i, K in enumerate(K_eq_list) if K > 100]
                                        if len(favorable_steps) > 1:
                                            st.write("**Additional driving forces:**")
                                            for idx in favorable_steps:
                                                if idx != max_K_idx:
                                                    st.write(f"• r{idx+1}: K = {K_eq_list[idx]:.2e}")
                                    
                                    # Practical recommendations
                                    st.subheader("💡 Practical Recommendations")
                                    
                                    recommendations = []
                                    
                                    # Overall thermodynamics recommendations
                                    if K_overall > 1000:
                                        recommendations.append("✅ **Excellent thermodynamics** - Reaction should proceed well under standard conditions")
                                        recommendations.append("🎯 **Focus on kinetics** - Thermodynamics are favorable, optimize rate constants and mass transfer")
                                    elif K_overall > 1:
                                        recommendations.append("✅ **Good thermodynamics** - Reaction is thermodynamically favorable")
                                        recommendations.append("⚖️ **Balanced optimization** - Consider both thermodynamic and kinetic improvements")
                                    elif K_overall > 0.001:
                                        recommendations.append("⚠️ **Marginal thermodynamics** - Consider temperature, pressure, or concentration optimization")
                                        recommendations.append("🌡️ **Temperature sensitivity** - Higher temperatures might improve equilibrium position")
                                    else:
                                        recommendations.append("❌ **Poor thermodynamics** - Reaction may require extreme conditions or process coupling")
                                        recommendations.append("🔄 **Process redesign** - Consider alternative pathways or reaction coupling")
                                    
                                    # Specific step recommendations
                                    unfavorable_steps = [i for i, K in enumerate(K_eq_list) if K < 0.01]
                                    if unfavorable_steps:
                                        step_names = [f"r{i+1}" for i in unfavorable_steps]
                                        recommendations.append(f"🔧 **Optimize specific steps:** {', '.join(step_names)} have very low K_eq values")
                                        recommendations.append("💡 **Consider:** Alternative catalysts, modified conditions, or step elimination")
                                    
                                    # Temperature optimization suggestions
                                    if any(K > 1000 for K in K_eq_list) and any(K < 0.001 for K in K_eq_list):
                                        recommendations.append("🌡️ **Temperature optimization** may help balance competing thermodynamic requirements")
                                        recommendations.append("📊 **Consider:** Van 't Hoff analysis to determine optimal temperature")
                                    
                                    # Process integration recommendations
                                    if 0.1 < K_overall < 10:
                                        recommendations.append("🔗 **Process integration** - Consider coupling with thermodynamically favorable reactions")
                                        recommendations.append("♻️ **Separation strategy** - Product removal might shift equilibrium favorably")
                                    
                                    # Display recommendations
                                    for i, rec in enumerate(recommendations, 1):
                                        st.write(f"{i}. {rec}")
                                    
                                    # Download section
                                    st.subheader("💾 Download Results")
                                    
                                    # Create comprehensive results for download
                                    download_data = []
                                    for i in range(len(K_eq_list)):
                                        if 'stoich' in loaded_files and i < len(loaded_files['stoich']):
                                            reaction_name = loaded_files['stoich'].iloc[i, 0]
                                        else:
                                            reaction_name = f"r{i+1}"
                                        
                                        download_data.append({
                                            'Reaction_ID': f"r{i+1}",
                                            'Reaction_Name': reaction_name,
                                            'k_forward': k_forward[i],
                                            'k_reverse': k_reverse[i],
                                            'K_equilibrium': K_eq_list[i],
                                            'Sigma': st.session_state.sigma_values[i],
                                            'K_eq_to_sigma': K_eq_list[i]**st.session_state.sigma_values[i],
                                            'ln_K_eq': math.log(K_eq_list[i]) if K_eq_list[i] > 0 else 'Undefined',
                                            'Delta_G_individual_kJ_mol': -R * T * math.log(K_eq_list[i]) / 1000 if K_eq_list[i] > 0 else 'Undefined',
                                            'Thermodynamic_Favorability': 'Forward' if K_eq_list[i] > 1 else 'Reverse' if K_eq_list[i] < 1 else 'Balanced'
                                        })
                                    
                                    # Add overall summary
                                    overall_summary = {
                                        'Reaction_ID': 'OVERALL',
                                        'Reaction_Name': 'Combined Process',
                                        'k_forward': '',
                                        'k_reverse': '',
                                        'K_equilibrium': K_overall,
                                        'Sigma': '',
                                        'K_eq_to_sigma': '',
                                        'ln_K_eq': math.log(K_overall) if K_overall > 0 else 'Undefined',
                                        'Delta_G_individual_kJ_mol': -R * T * math.log(K_overall) / 1000 if K_overall > 0 else 'Undefined',
                                        'Thermodynamic_Favorability': 'Forward' if K_overall > 1 else 'Reverse' if K_overall < 1 else 'Balanced'
                                    }
                                    download_data.append(overall_summary)
                                    
                                    # Add system parameters
                                    system_params = {
                                        'Reaction_ID': 'SYSTEM_PARAMS',
                                        'Reaction_Name': f'T={T}K, R={R}J/(mol·K)',
                                        'k_forward': f'Sum_sigma={sum(st.session_state.sigma_values):.3f}',
                                        'k_reverse': f'Num_reactions={len(K_eq_list)}',
                                        'K_equilibrium': '',
                                        'Sigma': '',
                                        'K_eq_to_sigma': '',
                                        'ln_K_eq': '',
                                        'Delta_G_individual_kJ_mol': '',
                                        'Thermodynamic_Favorability': ''
                                    }
                                    download_data.append(system_params)
                                    
                                    results_df = pd.DataFrame(download_data)
                                    csv = results_df.to_csv(index=False)
                                    
                                    st.download_button(
                                        label="📄 Download Equilibrium Analysis CSV",
                                        data=csv,
                                        file_name="uploaded_files_equilibrium_analysis.csv",
                                        mime="text/csv",
                                        help="Download comprehensive equilibrium analysis results with thermodynamic properties"
                                    )
                                    
                                except Exception as e:
                                    st.error(f"❌ Error in equilibrium analysis: {str(e)}")
                                    with st.expander("🔧 Error Details"):
                                        st.write(str(e))
                                        st.write("**Troubleshooting:**")
                                        st.write("• Verify all rate constants are positive numbers")
                                        st.write("• Check that forward and reverse rate constants are properly paired")
                                        st.write("• Ensure temperature is reasonable (200-2000 K)")
                                        st.write("• Validate stoichiometric numbers are appropriate")
                    
                    else:
                        st.warning("Stoichiometric matrix is required to determine the number of reactions.")
                
                else:
                    st.warning("Parameters file is required for thermodynamics analysis.")
