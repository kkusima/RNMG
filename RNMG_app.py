import streamlit as st
import pandas as pd
import numpy as np
import re
from io import StringIO
from collections import defaultdict, Counter
import math

# Page configuration
st.set_page_config(
    page_title="Reaction Network Matrix Generator",
    page_icon="🔗",
    layout="wide"
)

st.title("Reaction Network Matrix Generator")
st.markdown("Create stoichiometric and atomic matrices from chemical reaction mechanisms")

# Initialize session state
if 'reactions' not in st.session_state:
    st.session_state.reactions = []
if 'species_list' not in st.session_state:
    st.session_state.species_list = []
if 'atomic_matrix' not in st.session_state:
    st.session_state.atomic_matrix = pd.DataFrame()
if 'stoich_matrix' not in st.session_state:
    st.session_state.stoich_matrix = pd.DataFrame()
if 'param_data' not in st.session_state:
    st.session_state.param_data = []
if 'matrix_consistency' not in st.session_state:
    st.session_state.matrix_consistency = {'consistent': True, 'atomic_species': [], 'stoich_species': []}
if 'mass_balance_result' not in st.session_state:
    st.session_state.mass_balance_result = {'balanced': None, 'message': ''}
if 'app_mode' not in st.session_state:
    st.session_state.app_mode = 'Generator'
if 'uploaded_files' not in st.session_state:
    st.session_state.uploaded_files = {'atomic': None, 'stoich': None, 'params': None}
if 'analysis_mode' not in st.session_state:
    st.session_state.analysis_mode = 'Mass Balance'
if 'sigma_values' not in st.session_state:
    st.session_state.sigma_values = []
if 'equilibrium_constants' not in st.session_state:
    st.session_state.equilibrium_constants = []

def parse_species(species_str):
    """Parse a species string to extract species name and stoichiometric coefficient"""
    species_str = species_str.strip()
    
    # Handle coefficients like 2O*, 2*
    coeff_match = re.match(r'^(\d+)(.+)$', species_str)
    if coeff_match:
        coeff = int(coeff_match.group(1))
        species = coeff_match.group(2)
    else:
        coeff = 1
        species = species_str
    
    return species, coeff

def parse_reaction(reaction_str):
    """Parse a reaction string and extract reactants and products"""
    # Check if the required ⇌ arrow is present
    if '⇌' not in reaction_str:
        st.error(f"❌ Please use only ⇌ arrows. Copy the arrow from the box above. Your input: {reaction_str}")
        return None, None, None
    
    # Split into reactants and products
    parts = reaction_str.split('⇌')
    if len(parts) != 2:
        st.error(f"Invalid reaction format: {reaction_str}")
        return None, None, None
    
    reactants_str, products_str = parts
    
    # Parse reactants
    reactants = {}
    for reactant in reactants_str.split('+'):
        reactant = reactant.strip()
        if reactant:
            species, coeff = parse_species(reactant)
            reactants[species] = reactants.get(species, 0) + coeff
    
    # Parse products
    products = {}
    for product in products_str.split('+'):
        product = product.strip()
        if product:
            species, coeff = parse_species(product)
            products[species] = products.get(species, 0) + coeff
    
    return reactants, products, reaction_str

def get_atomic_composition():
    """Define atomic composition for common species"""
    compositions = {
        'CO': {'C': 1, 'O': 1},
        'O2': {'O': 2},
        'CO2': {'C': 1, 'O': 2},
        'CO*': {'C': 1, 'O': 1, '*': 1},
        'O*': {'O': 1, '*': 1},
        'O2*': {'O': 2, '*': 1},
        '*': {'*': 1},
        'H2': {'H': 2},
        'H2O': {'H': 2, 'O': 1},
        'CH4': {'C': 1, 'H': 4},
        'NH3': {'N': 1, 'H': 3},
        'N2': {'N': 2},
        'H*': {'H': 1, '*': 1},
        'OH*': {'O': 1, 'H': 1, '*': 1},
        'C*': {'C': 1, '*': 1},
        'N*': {'N': 1, '*': 1}
    }
    return compositions

def get_species_atomic_composition(species):
    """Get atomic composition for a species, handling gas phase notation"""
    compositions = get_atomic_composition()
    
    # Remove gas phase notation (g), (l), (s) if present
    clean_species = re.sub(r'\([gls]\)', '', species)
    
    # Look up composition
    if clean_species in compositions:
        return compositions[clean_species]
    else:
        # If not found, try to parse simple species
        # This is a basic parser - can be expanded for more complex species
        return parse_species_formula(clean_species)

def parse_species_formula(formula):
    """Basic parser for simple chemical formulas"""
    composition = {}
    
    # Handle surface sites
    if formula == '*':
        return {'*': 1}
    
    # Handle adsorbed species (ending with *)
    if formula.endswith('*'):
        base_formula = formula[:-1]
        composition['*'] = 1
        # Parse the base formula
        base_comp = parse_species_formula(base_formula)
        composition.update(base_comp)
        return composition
    
    # Simple regex pattern to match element-number pairs
    pattern = r'([A-Z][a-z]?)(\d*)'
    matches = re.findall(pattern, formula)
    
    for element, count in matches:
        count = int(count) if count else 1
        composition[element] = composition.get(element, 0) + count
    
    return composition if composition else {'Unknown': 1}

def update_matrices():
    """Update atomic and stoichiometric matrices based on current reactions"""
    if not st.session_state.reactions:
        st.session_state.atomic_matrix = pd.DataFrame()
        st.session_state.stoich_matrix = pd.DataFrame()
        st.session_state.species_list = []
        return
    
    # Collect all species
    all_species = set()
    for _, reactants, products, _ in st.session_state.reactions:
        all_species.update(reactants.keys())
        all_species.update(products.keys())
    
    # Separate and sort gas and surface species, with empty sites (*) last
    gas_species = sorted([s for s in all_species if not s.endswith('*')])
    surface_species = [s for s in all_species if s.endswith('*')]
    
    # Sort surface species with empty sites (*) at the end
    other_surface = sorted([s for s in surface_species if s != '*'])
    empty_sites = [s for s in surface_species if s == '*']
    surface_species = other_surface + empty_sites
    
    # Create consistent species order: gas species first, then surface species, with * last
    ordered_species = gas_species + surface_species
    st.session_state.species_list = ordered_species
    
    # Create atomic matrix with consistent ordering
    all_atoms = set()
    
    # Collect all atom types using the new function
    for species in ordered_species:
        composition = get_species_atomic_composition(species)
        all_atoms.update(composition.keys())
    
    # Sort atoms with "*" (surface sites) at the end
    other_atoms = sorted([atom for atom in all_atoms if atom != '*'])
    surface_atom = [atom for atom in all_atoms if atom == '*']
    all_atoms = other_atoms + surface_atom
    
    # Build atomic matrix with consistent species ordering
    atomic_data = []
    for atom in all_atoms:
        row = [atom]
        for species in ordered_species:  # Use ordered_species to ensure consistent ordering
            composition = get_species_atomic_composition(species)
            row.append(composition.get(atom, 0))
        atomic_data.append(row)
    
    columns = ['A\\S'] + ordered_species
    st.session_state.atomic_matrix = pd.DataFrame(atomic_data, columns=columns)
    
    # Create stoichiometric matrix with same species ordering
    stoich_data = []
    
    # Get gas and surface species in the same order as atomic matrix
    gas_species = [s for s in ordered_species if not s.endswith('*')]
    surface_species = [s for s in ordered_species if s.endswith('*')]
    
    # Column headers: P_ for gas phase, theta_ for surface species (with * last)
    # Use the same ordering as atomic matrix
    stoich_columns = ['r\\S'] + [f'P_{s}' for s in gas_species] + [f'theta_{s}' for s in surface_species]
    
    for i, (_, reactants, products, _) in enumerate(st.session_state.reactions):
        row = [f'r{i+1}']
        
        # Gas phase species (in same order as atomic matrix)
        for species in gas_species:
            net_coeff = products.get(species, 0) - reactants.get(species, 0)
            row.append(net_coeff)
        
        # Surface species (in same order as atomic matrix)
        for species in surface_species:
            net_coeff = products.get(species, 0) - reactants.get(species, 0)
            row.append(net_coeff)
        
        stoich_data.append(row)
    
    st.session_state.stoich_matrix = pd.DataFrame(stoich_data, columns=stoich_columns)
    
    # Verification check: ensure both matrices have consistent species ordering
    verify_matrix_consistency()

def verify_matrix_consistency():
    """Verify that atomic and stoichiometric matrices have consistent species ordering"""
    if st.session_state.atomic_matrix.empty or st.session_state.stoich_matrix.empty:
        return True
    
    # Get species from atomic matrix (excluding first column A\S)
    atomic_species = st.session_state.atomic_matrix.columns[1:].tolist()
    
    # Get species from stoichiometric matrix (excluding first column r\S)
    # Remove prefixes P_ and theta_ to get actual species names
    stoich_species = []
    for col in st.session_state.stoich_matrix.columns[1:]:
        if col.startswith('P_'):
            stoich_species.append(col[2:])  # Remove 'P_' prefix
        elif col.startswith('theta_'):
            stoich_species.append(col[6:])  # Remove 'theta_' prefix
    
    # Check if species lists match
    consistency_check = atomic_species == stoich_species
    
    # Store consistency status in session state for display
    st.session_state.matrix_consistency = {
        'consistent': consistency_check,
        'atomic_species': atomic_species,
        'stoich_species': stoich_species
    }
    
    return consistency_check

def check_mass_balance(atomic_matrix=None, stoich_matrix=None):
    """Check if mass is conserved across all reactions using A*v=0"""
    if atomic_matrix is None:
        atomic_matrix = st.session_state.atomic_matrix
    if stoich_matrix is None:
        stoich_matrix = st.session_state.stoich_matrix
        
    if atomic_matrix.empty or stoich_matrix.empty:
        return False, "Matrices not available for mass balance check"
    
    # Extract atomic matrix (excluding first column A\S)
    atomic_matrix_vals = atomic_matrix.iloc[:, 1:].values
    
    # Extract stoichiometric matrix (excluding first column r\S)
    stoich_matrix_vals = stoich_matrix.iloc[:, 1:].values
    
    errors = []
    mass_balanced = True
    
    # Check mass balance for each reaction
    for i, reaction_row in enumerate(stoich_matrix_vals):
        # Compute A * v for reaction i
        mass_balance_result = np.dot(atomic_matrix_vals, reaction_row)
        
        # Check if all elements are approximately zero (mass balanced)
        if not np.allclose(mass_balance_result, 0, atol=1e-10):
            mass_balanced = False
            reaction_name = stoich_matrix.iloc[i, 0]
            errors.append({
                'reaction': i + 1,
                'reaction_name': reaction_name,
                'imbalance': mass_balance_result.tolist()
            })
    
    if mass_balanced:
        return True, "✅ Mass is conserved across all reactions."
    else:
        error_message = "❌ Mass is NOT conserved in the following reactions:\n"
        for error in errors:
            error_message += f"• Reaction {error['reaction']} ({error['reaction_name']}): Imbalance = {error['imbalance']}\n"
        error_message += "\nPlease check and correct the Atomic and/or Stoichiometric matrices."
        return False, error_message

def create_default_parameters():
    """Create default parameter entries for the current reaction set"""
    params = []
    
    # Global parameters
    params.append({'Reaction_Descrp': '', 'Parameter': 'T', 'Values': 320.0, 'Units': 'K'})
    params.append({'Reaction_Descrp': '', 'Parameter': 'R', 'Values': 8.31446, 'Units': 'JK^-1mol^-1'})
    
    # Pressure parameters for gas species (in consistent order)
    gas_species = [s for s in st.session_state.species_list if not s.endswith('*')]
    for i, species in enumerate(gas_species):
        params.append({'Reaction_Descrp': species, 'Parameter': f'P{i+1}', 'Values': 1.0e-8, 'Units': 'bar'})
    
    # Rate constants for reactions
    for i in range(len(st.session_state.reactions)):
        reaction_id = f'r{i+1}'
        params.append({'Reaction_Descrp': reaction_id, 'Parameter': f'k{i+1}f', 'Values': 1.0, 'Units': '-'})
        params.append({'Reaction_Descrp': '', 'Parameter': f'k{i+1}r', 'Values': 1.0, 'Units': '-'})
    
    # Constants
    for i in range(len(st.session_state.reactions)):
        for direction in ['f', 'r']:
            for const_type in ['a', 'b', 'c']:
                params.append({'Reaction_Descrp': 'const', 'Parameter': f'{const_type}{i+1}{direction}', 'Values': 1.0, 'Units': '-'})
    
    return params

def load_uploaded_files():
    """Load and process uploaded files"""
    files_loaded = {}
    
    if st.session_state.uploaded_files['atomic'] is not None:
        try:
            atomic_df = pd.read_csv(st.session_state.uploaded_files['atomic'])
            files_loaded['atomic'] = atomic_df
        except Exception as e:
            st.error(f"Error loading atomic matrix: {e}")
            return None
    
    if st.session_state.uploaded_files['stoich'] is not None:
        try:
            stoich_df = pd.read_csv(st.session_state.uploaded_files['stoich'])
            files_loaded['stoich'] = stoich_df
        except Exception as e:
            st.error(f"Error loading stoichiometric matrix: {e}")
            return None
    
    if st.session_state.uploaded_files['params'] is not None:
        try:
            params_df = pd.read_csv(st.session_state.uploaded_files['params'])
            files_loaded['params'] = params_df
        except Exception as e:
            st.error(f"Error loading parameters: {e}")
            return None
    
    return files_loaded

def calculate_equilibrium_constants(params_df, sigma_values, k_pattern="auto"):
    """Calculate individual and overall equilibrium constants"""
    try:
        # Extract temperature and R
        T_row = params_df[params_df['Parameter'] == 'T']
        R_row = params_df[params_df['Parameter'] == 'R']
        
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
        
        # Find all rate constants - be flexible with naming
        k_forward = []
        k_reverse = []
        
        # Get all parameters that start with 'k' and contain numbers
        k_params = params_df[params_df['Parameter'].str.contains(r'^k\d+', regex=True, na=False)]
        k_params = k_params.sort_values('Parameter')
        
        if k_pattern == "alternating":
            # Every other k is reverse (k1=forward, k2=reverse, k3=forward, k4=reverse...)
            k_values = k_params['Values'].tolist()
            for i in range(0, len(k_values), 2):
                if i < len(k_values):
                    k_forward.append(k_values[i])
                if i + 1 < len(k_values):
                    k_reverse.append(k_values[i + 1])
                else:
                    k_reverse.append(1.0)  # Default if missing
                    
        else:
            # Try to detect pattern automatically
            # Look for f/r suffixes first
            for i in range(len(sigma_values)):
                # Try different naming conventions
                possible_forward = [f'k{i+1}f', f'k{i+1}F', f'k_{i+1}_f', f'k{i+1}', f'k{2*i+1}']
                possible_reverse = [f'k{i+1}r', f'k{i+1}R', f'k_{i+1}_r', f'k{i+1}_rev', f'k{2*i+2}']
                
                kf_found = False
                kr_found = False
                
                for kf_name in possible_forward:
                    kf_row = params_df[params_df['Parameter'] == kf_name]
                    if not kf_row.empty:
                        k_forward.append(kf_row['Values'].iloc[0])
                        kf_found = True
                        break
                
                for kr_name in possible_reverse:
                    kr_row = params_df[params_df['Parameter'] == kr_name]
                    if not kr_row.empty:
                        k_reverse.append(kr_row['Values'].iloc[0])
                        kr_found = True
                        break
                
                if not kf_found:
                    st.warning(f"Forward rate constant for reaction {i+1} not found, using 1.0")
                    k_forward.append(1.0)
                    
                if not kr_found:
                    st.warning(f"Reverse rate constant for reaction {i+1} not found, using 1.0")
                    k_reverse.append(1.0)
        
        # Ensure we have the right number of rate constants
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
            if K_i > 0:  # Ensure positive value for logarithm
                K_overall *= K_i ** sigma_i
            else:
                st.warning(f"Non-positive equilibrium constant for reaction {i+1}")
                K_overall *= 1.0
        
        return K_eq, K_overall, T, R, k_forward, k_reverse
        
    except Exception as e:
        st.error(f"Error calculating equilibrium constants: {e}")
        return [], 1.0, 300.0, 8.314, [], []

# Mode selector in sidebar
with st.sidebar:
    st.header("Application Mode")
    mode = st.radio(
        "Select Mode:",
        ["Generator", "Analysis"],
        index=0 if st.session_state.app_mode == 'Generator' else 1
    )
    st.session_state.app_mode = mode

# Generator Mode (Original functionality)
if st.session_state.app_mode == 'Generator':
    # Original sidebar content for reaction input
    with st.sidebar:
        st.header("Reaction Input")
        
        # Reminder about arrow and single reaction
        st.markdown("**Use ⇌ arrow (copy above) • One reaction at a time**")
        
        reaction_input = st.text_area(
            "Enter ONE reaction:",
            placeholder="CO(g) + * ⇌ CO*",
            help="Use ⇌ for arrows. Use * for surface sites. Add ONE reaction, then click 'Add Reaction'.",
            height=80
        )
        
        if st.button("Add Reaction", type="primary"):
            if reaction_input.strip():
                reactants, products, original = parse_reaction(reaction_input.strip())
                if reactants is not None and products is not None:
                    st.session_state.reactions.append((len(st.session_state.reactions), reactants, products, original))
                    update_matrices()
                    st.session_state.param_data = create_default_parameters()
                    st.session_state.mass_balance_result = {'balanced': None, 'message': ''}  # Reset mass balance
                    st.success("Reaction added!")
                    st.rerun()
        
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
                        st.session_state.mass_balance_result = {'balanced': None, 'message': ''}  # Reset mass balance
                        st.rerun()
        
        if st.button("Clear All", type="secondary"):
            st.session_state.reactions = []
            st.session_state.species_list = []
            st.session_state.atomic_matrix = pd.DataFrame()
            st.session_state.stoich_matrix = pd.DataFrame()
            st.session_state.param_data = []
            st.session_state.mass_balance_result = {'balanced': None, 'message': ''}
            st.rerun()

    # Main content area for Generator mode
    if st.session_state.reactions:
        # Tabs for different matrices
        tab1, tab2, tab3, tab4 = st.tabs(["📊 Stoichiometric Matrix", "⚛️ Atomic Matrix", "⚙️ Parameters", "🌡️ Thermodynamics"])
        
        with tab1:
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
        
        with tab2:
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
        
        with tab3:
            st.subheader("Parameter Table")
            
            if st.session_state.param_data:
                # Create editable dataframe
                param_df = pd.DataFrame(st.session_state.param_data)
                
                # Display editable table
                edited_df = st.data_editor(
                    param_df,
                    use_container_width=True,
                    num_rows="dynamic",
                    column_config={
                        "Reaction_Descrp": st.column_config.TextColumn("Reaction Description"),
                        "Parameter": st.column_config.TextColumn("Parameter"),
                        "Values": st.column_config.NumberColumn("Values", format="%.2e"),
                        "Units": st.column_config.TextColumn("Units")
                    }
                )
                
                # Update session state with edits
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
        
        with tab4:
            st.subheader("Thermodynamics Analysis")
            
            if st.session_state.param_data:
                # Convert param_data to DataFrame for processing
                params_df = pd.DataFrame(st.session_state.param_data)
                
                # Parameter inspection section
                with st.expander("🔍 Inspect Generated Parameters"):
                    st.write("**Current Parameters:**")
                    st.dataframe(params_df, use_container_width=True)
                    
                    # Highlight rate constants
                    k_params = params_df[params_df['Parameter'].str.contains(r'^k\d+', regex=True, na=False)]
                    if not k_params.empty:
                        st.write("**Rate Constants Found:**")
                        st.dataframe(k_params, use_container_width=True)
                
                # Rate constant pattern selector
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
                
                # Edit parameters inline
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
                
                # Update session state with edits
                st.session_state.param_data = edited_params_gen.to_dict('records')
                
                # Number of reactions
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
                
                # Calculate equilibrium constants
                if st.button("🧮 Calculate Overall Equilibrium", type="primary", key="gen_calc_eq"):
                    K_eq_list_gen, K_overall_gen, T_gen, R_gen, k_forward_gen, k_reverse_gen = calculate_equilibrium_constants(
                        edited_params_gen, st.session_state.sigma_values_gen, k_pattern_gen
                    )
                    
                    st.markdown("---")
                    st.subheader("🎯 Results")
                    
                    # Display rate constants used
                    st.markdown("**Rate Constants Used:**")
                    rate_data_gen = []
                    for i, (kf, kr) in enumerate(zip(k_forward_gen, k_reverse_gen)):
                        # Get original reaction for reference
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
                    
                    # Display equilibrium calculation
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
                    
                    # Display overall equilibrium constant
                    st.markdown("**Overall Equilibrium Constant:**")
                    st.latex(r"K_{overall,eq} = \prod_{i} K_{i,eq}^{\sigma_i}")
                    
                    # Format scientific notation for display
                    if K_overall_gen != 0:
                        exponent = math.floor(math.log10(abs(K_overall_gen)))
                        mantissa = K_overall_gen / (10 ** exponent)
                        st.markdown(f"**K_overall = {mantissa:.4f} × 10^{exponent}**")
                    else:
                        st.markdown("**K_overall = 0**")
                    
                    # Additional thermodynamic info
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
    
        # Matrix consistency check
        st.markdown("---")
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
        
        # Mass Balance Check Section
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
        
        # Summary section
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
        # Welcome message and example
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

# Analysis Mode (New functionality)
else:
    with st.sidebar:
        st.header("File Upload")
        
        # File uploaders
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
        
        # Analysis mode selector
        st.header("Analysis Type")
        analysis_mode = st.radio(
            "Select Analysis:",
            ["Mass Balance", "Thermodynamics"],
            index=0 if st.session_state.analysis_mode == 'Mass Balance' else 1
        )
        st.session_state.analysis_mode = analysis_mode

    # Main content for Analysis mode
    st.markdown("## 🔬 Analysis Mode")
    
    # Check if files are uploaded
    files_uploaded = any(st.session_state.uploaded_files[key] is not None for key in st.session_state.uploaded_files)
    
    if not files_uploaded:
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
        # Load files
        loaded_files = load_uploaded_files()
        
        if loaded_files:
            # Display loaded files info
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
                    # Display matrices
                    col1, col2 = st.columns(2)
                    
                    with col1:
                        st.subheader("Atomic Matrix")
                        st.dataframe(loaded_files['atomic'], use_container_width=True)
                    
                    with col2:
                        st.subheader("Stoichiometric Matrix")
                        st.dataframe(loaded_files['stoich'], use_container_width=True)
                    
                    # Mass balance check
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
                    # Display and edit parameters
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
                            
                            # Check for different patterns
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
                    
                    # Make parameters editable
                    edited_params = st.data_editor(
                        loaded_files['params'],
                        use_container_width=True,
                        num_rows="dynamic",
                        column_config={
                            "Values": st.column_config.NumberColumn("Values", format="%.6e")
                        }
                    )
                    
                    # Determine number of reactions
                    if 'stoich' in loaded_files:
                        num_reactions = loaded_files['stoich'].shape[0]
                        
                        st.markdown("---")
                        st.subheader("Stoichiometric Numbers (σ)")
                        st.markdown("Enter the stoichiometric number for each reaction step:")
                        
                        # Initialize sigma values if not set
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
                            
                            # Display individual equilibrium constants
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
                            
                            # Format scientific notation for display
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
