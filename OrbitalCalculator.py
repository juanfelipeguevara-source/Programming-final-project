import tkinter as tk
from tkinter import ttk, messagebox
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

# --- Constantes ---
R_EARTH_KM = 6378.137  # Radio terrestre en km
MU_EARTH = 398600.4418  # Parámetro gravitacional terrestre en km³/s²
MU_SUN = 1.32712440018e11  # km³/s²
MU_MOON = 4902.8  # km³/s²
J2 = 1.08263e-3  # Coeficiente de aplanamiento J2
EARTH_ROT = 7.2921159e-5  # Velocidad de rotación terrestre (rad/s)
AU = 149597870.7  # Unidad Astronómica en km
MOON_DIST = 384400  # Distancia promedio Tierra-Luna en km
GEO_ALT = 35786  # Altitud órbita geoestacionaria en km
MOON_SOI = 66100  # Esfera de influencia lunar aproximada en km

# --- Funciones auxiliares ---
def calc_period(a_km):
    """Calcula el periodo orbital en segundos"""
    return 2 * np.pi * np.sqrt(a_km**3 / MU_EARTH)

def calc_velocity(r_km, a_km):
    """Calcula la velocidad en una órbita en km/s"""
    return np.sqrt(MU_EARTH * (2/r_km - 1/a_km))

def calc_perturbation_j2(a, e, i):
    """Calcula la precesión nodal y del perigeo debido a J2"""
    n = np.sqrt(MU_EARTH / a**3)  # Movimiento medio
    p = a * (1 - e**2)  # Parámetro orbital
    
    # Precesión del nodo ascendente (rad/s)
    dOmega_dt = -1.5 * n * J2 * (R_EARTH_KM / p)**2 * np.cos(i)
    
    # Precesión del argumento del perigeo (rad/s)
    domega_dt = 1.5 * n * J2 * (R_EARTH_KM / p)**2 * (2 - 2.5 * np.sin(i)**2)
    
    return dOmega_dt, domega_dt

def calc_third_body_perturbation(r_sat, r_body, mu_body):
    """Calcula la aceleración debido a perturbación de tercer cuerpo"""
    r_sat_body = r_body - r_sat
    r_sb_norm = np.linalg.norm(r_sat_body)
    r_b_norm = np.linalg.norm(r_body)
    
    accel = mu_body * (r_sat_body / r_sb_norm**3 - r_body / r_b_norm**3)
    return accel

def atmospheric_density(h_km):
    """Modelo exponencial simplificado de densidad atmosférica"""
    if h_km > 1000:
        return 0
    elif h_km > 500:
        rho0 = 6e-16
        H = 60
        h0 = 500
    elif h_km > 300:
        rho0 = 1.5e-13
        H = 50
        h0 = 300
    elif h_km > 200:
        rho0 = 2.5e-11
        H = 40
        h0 = 200
    else:
        rho0 = 5e-10
        H = 30
        h0 = 150
    
    return rho0 * np.exp(-(h_km - h0) / H)

def hohmann_transfer(r1, r2):
    """Calcula una transferencia de Hohmann"""
    a_initial = r1
    a_target = r2
    a_transfer = (r1 + r2) / 2
    
    v1_initial = calc_velocity(r1, a_initial)
    v1_transfer = calc_velocity(r1, a_transfer)
    dv1 = abs(v1_transfer - v1_initial)
    
    v2_transfer = calc_velocity(r2, a_transfer)
    v2_circular = calc_velocity(r2, a_target)
    dv2 = abs(v2_circular - v2_transfer)
    
    T_transfer = calc_period(a_transfer) / 2
    
    return dv1, dv2, dv1 + dv2, T_transfer

def bielliptic_transfer(r1, r2, r_intermediate):
    """Calcula una transferencia bi-elíptica"""
    # Primera transferencia: r1 -> r_intermediate
    a_transfer1 = (r1 + r_intermediate) / 2
    v1_initial = calc_velocity(r1, r1)
    v1_transfer1 = calc_velocity(r1, a_transfer1)
    dv1 = abs(v1_transfer1 - v1_initial)
    
    # Segunda maniobra en r_intermediate
    v_inter_transfer1 = calc_velocity(r_intermediate, a_transfer1)
    a_transfer2 = (r_intermediate + r2) / 2
    v_inter_transfer2 = calc_velocity(r_intermediate, a_transfer2)
    dv2 = abs(v_inter_transfer2 - v_inter_transfer1)
    
    # Tercera maniobra en r2
    v2_transfer2 = calc_velocity(r2, a_transfer2)
    v2_circular = calc_velocity(r2, r2)
    dv3 = abs(v2_circular - v2_transfer2)
    
    T_transfer1 = calc_period(a_transfer1) / 2
    T_transfer2 = calc_period(a_transfer2) / 2
    T_total = T_transfer1 + T_transfer2
    
    return dv1, dv2, dv3, dv1 + dv2 + dv3, T_total

def calc_plane_change_dv(v, delta_i):
    """Calcula el delta-v para cambio de plano orbital"""
    return 2 * v * np.sin(delta_i / 2)

def get_orbit_classification(altitude_km):
    """Clasifica el tipo de órbita según su altitud"""
    if altitude_km < 2000:
        return "LEO (Low Earth Orbit)"
    elif altitude_km < 35000:
        return "MEO (Medium Earth Orbit)"
    elif 35000 <= altitude_km <= 36500:
        return "GEO (Geostationary Orbit)"
    elif altitude_km < 100000:
        return "HEO (High Earth Orbit)"
    elif altitude_km < MOON_DIST * 0.5:
        return "Cislunar (Tierra-Luna)"
    else:
        return "Órbita Lunar/Translunar"

class OrbitalCalculatorApp:
    def __init__(self, root):
        self.root = root
        self.root.title("🛰️ Calculadora Orbital Avanzada - Órbitas Extendidas")
        self.root.geometry("1600x900")
        
        # Variables con rangos extendidos
        self.alt_p_inicial = tk.DoubleVar(value=400.0)
        self.alt_a_inicial = tk.DoubleVar(value=600.0)
        self.alt_target = tk.DoubleVar(value=35786.0)  # GEO por defecto
        self.sat_mass = tk.DoubleVar(value=100.0)
        self.inclination_initial = tk.DoubleVar(value=28.5)
        self.inclination_target = tk.DoubleVar(value=0.0)  # GEO ecuatorial
        self.maneuver_type = tk.StringVar(value="Hohmann")
        self.alt_intermediate = tk.DoubleVar(value=50000.0)
        self.include_j2 = tk.BooleanVar(value=True)
        self.include_drag = tk.BooleanVar(value=False)
        self.include_sun = tk.BooleanVar(value=True)  # Importante para órbitas altas
        self.include_moon = tk.BooleanVar(value=True)  # Importante para órbitas altas
        self.drag_coef = tk.DoubleVar(value=2.2)
        self.cross_section = tk.DoubleVar(value=1.0)
        self.orbit_preset = tk.StringVar(value="Custom")
        
        # Variables para almacenar resultados
        self.current_results = None
        
        self.create_widgets()
        
    def create_widgets(self):
        # Frame principal
        main_frame = ttk.Frame(self.root)
        main_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        # Panel izquierdo (Inputs)
        left_frame = ttk.LabelFrame(main_frame, text="⚙️ Parámetros de Entrada", padding=10)
        left_frame.grid(row=0, column=0, sticky="nsew", padx=(0, 5))
        
        # Scroll para inputs
        canvas = tk.Canvas(left_frame, width=350)
        scrollbar = ttk.Scrollbar(left_frame, orient="vertical", command=canvas.yview)
        scrollable_frame = ttk.Frame(canvas)
        
        scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
        )
        
        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)
        
        row = 0
        
        # Órbitas Preconfiguradas
        ttk.Label(scrollable_frame, text="🎯 ÓRBITAS PRECONFIGURADAS", font=('Arial', 10, 'bold')).grid(row=row, column=0, columnspan=3, pady=5)
        row += 1
        
        preset_frame = ttk.Frame(scrollable_frame)
        preset_frame.grid(row=row, column=0, columnspan=3, pady=5)
        
        presets = [
            ("LEO (400km)", "LEO"),
            ("MEO/GPS (20200km)", "MEO"),
            ("GTO", "GTO"),
            ("GEO (35786km)", "GEO"),
            ("Lunar", "LUNAR"),
            ("Custom", "Custom")
        ]
        
        for i, (text, value) in enumerate(presets):
            ttk.Radiobutton(preset_frame, text=text, variable=self.orbit_preset, 
                          value=value, command=self.apply_preset).grid(row=i//2, column=i%2, sticky="w", padx=5)
        row += 1
        
        ttk.Separator(scrollable_frame, orient='horizontal').grid(row=row, column=0, columnspan=3, sticky='ew', pady=10)
        row += 1
        
        # Órbita Inicial
        ttk.Label(scrollable_frame, text="🌍 ÓRBITA INICIAL", font=('Arial', 10, 'bold')).grid(row=row, column=0, columnspan=3, pady=5)
        row += 1
        
        ttk.Label(scrollable_frame, text="Altitud Perigeo (km):").grid(row=row, column=0, sticky="w", pady=2)
        self.perigee_scale = ttk.Scale(scrollable_frame, from_=150, to=100000, variable=self.alt_p_inicial, orient=tk.HORIZONTAL, length=150)
        self.perigee_scale.grid(row=row, column=1, pady=2)
        perigee_entry = ttk.Entry(scrollable_frame, textvariable=self.alt_p_inicial, width=10)
        perigee_entry.grid(row=row, column=2, padx=2)
        row += 1
        
        ttk.Label(scrollable_frame, text="Altitud Apogeo (km):").grid(row=row, column=0, sticky="w", pady=2)
        self.apogee_scale = ttk.Scale(scrollable_frame, from_=150, to=400000, variable=self.alt_a_inicial, orient=tk.HORIZONTAL, length=150)
        self.apogee_scale.grid(row=row, column=1, pady=2)
        apogee_entry = ttk.Entry(scrollable_frame, textvariable=self.alt_a_inicial, width=10)
        apogee_entry.grid(row=row, column=2, padx=2)
        row += 1
        
        ttk.Label(scrollable_frame, text="Inclinación (°):").grid(row=row, column=0, sticky="w", pady=2)
        ttk.Scale(scrollable_frame, from_=0, to=90, variable=self.inclination_initial, orient=tk.HORIZONTAL, length=150).grid(row=row, column=1, pady=2)
        inc_initial_entry = ttk.Entry(scrollable_frame, textvariable=self.inclination_initial, width=10)
        inc_initial_entry.grid(row=row, column=2, padx=2)
        row += 1
        
        # Órbita Objetivo
        ttk.Separator(scrollable_frame, orient='horizontal').grid(row=row, column=0, columnspan=3, sticky='ew', pady=10)
        row += 1
        ttk.Label(scrollable_frame, text="🎯 ÓRBITA OBJETIVO", font=('Arial', 10, 'bold')).grid(row=row, column=0, columnspan=3, pady=5)
        row += 1
        
        ttk.Label(scrollable_frame, text="Altitud Circular (km):").grid(row=row, column=0, sticky="w", pady=2)
        self.target_scale = ttk.Scale(scrollable_frame, from_=500, to=400000, variable=self.alt_target, orient=tk.HORIZONTAL, length=150)
        self.target_scale.grid(row=row, column=1, pady=2)
        target_entry = ttk.Entry(scrollable_frame, textvariable=self.alt_target, width=10)
        target_entry.grid(row=row, column=2, padx=2)
        row += 1
        
        ttk.Label(scrollable_frame, text="Inclinación (°):").grid(row=row, column=0, sticky="w", pady=2)
        ttk.Scale(scrollable_frame, from_=0, to=90, variable=self.inclination_target, orient=tk.HORIZONTAL, length=150).grid(row=row, column=1, pady=2)
        inc_target_entry = ttk.Entry(scrollable_frame, textvariable=self.inclination_target, width=10)
        inc_target_entry.grid(row=row, column=2, padx=2)
        row += 1
        
        # Tipo de Maniobra
        ttk.Separator(scrollable_frame, orient='horizontal').grid(row=row, column=0, columnspan=3, sticky='ew', pady=10)
        row += 1
        ttk.Label(scrollable_frame, text="🚀 TIPO DE MANIOBRA", font=('Arial', 10, 'bold')).grid(row=row, column=0, columnspan=3, pady=5)
        row += 1
        
        ttk.Radiobutton(scrollable_frame, text="Hohmann (2 impulsos)", variable=self.maneuver_type, value="Hohmann").grid(row=row, column=0, columnspan=3, sticky="w")
        row += 1
        ttk.Radiobutton(scrollable_frame, text="Bi-elíptica (3 impulsos)", variable=self.maneuver_type, value="Bielliptic").grid(row=row, column=0, columnspan=3, sticky="w")
        row += 1
        
        ttk.Label(scrollable_frame, text="Altitud Intermedia (km):").grid(row=row, column=0, sticky="w", pady=2)
        ttk.Scale(scrollable_frame, from_=10000, to=500000, variable=self.alt_intermediate, orient=tk.HORIZONTAL, length=150).grid(row=row, column=1, pady=2)
        intermediate_entry = ttk.Entry(scrollable_frame, textvariable=self.alt_intermediate, width=10)
        intermediate_entry.grid(row=row, column=2, padx=2)
        row += 1
        
        # Perturbaciones
        ttk.Separator(scrollable_frame, orient='horizontal').grid(row=row, column=0, columnspan=3, sticky='ew', pady=10)
        row += 1
        ttk.Label(scrollable_frame, text="🌊 PERTURBACIONES", font=('Arial', 10, 'bold')).grid(row=row, column=0, columnspan=3, pady=5)
        row += 1
        
        ttk.Checkbutton(scrollable_frame, text="J2 (aplanamiento)", variable=self.include_j2).grid(row=row, column=0, columnspan=3, sticky="w")
        row += 1
        ttk.Checkbutton(scrollable_frame, text="Arrastre atmosférico", variable=self.include_drag).grid(row=row, column=0, columnspan=3, sticky="w")
        row += 1
        ttk.Checkbutton(scrollable_frame, text="Perturbación solar", variable=self.include_sun).grid(row=row, column=0, columnspan=3, sticky="w")
        row += 1
        ttk.Checkbutton(scrollable_frame, text="Perturbación lunar", variable=self.include_moon).grid(row=row, column=0, columnspan=3, sticky="w")
        row += 1
        
        ttk.Label(scrollable_frame, text="Coef. Arrastre (Cd):").grid(row=row, column=0, sticky="w", pady=2)
        ttk.Entry(scrollable_frame, textvariable=self.drag_coef, width=10).grid(row=row, column=1, sticky="w")
        row += 1
        
        ttk.Label(scrollable_frame, text="Sección Transv. (m²):").grid(row=row, column=0, sticky="w", pady=2)
        ttk.Entry(scrollable_frame, textvariable=self.cross_section, width=10).grid(row=row, column=1, sticky="w")
        row += 1
        
        # Satélite
        ttk.Separator(scrollable_frame, orient='horizontal').grid(row=row, column=0, columnspan=3, sticky='ew', pady=10)
        row += 1
        ttk.Label(scrollable_frame, text="🛰️ SATÉLITE", font=('Arial', 10, 'bold')).grid(row=row, column=0, columnspan=3, pady=5)
        row += 1
        
        ttk.Label(scrollable_frame, text="Masa (kg):").grid(row=row, column=0, sticky="w", pady=2)
        ttk.Entry(scrollable_frame, textvariable=self.sat_mass, width=10).grid(row=row, column=1, sticky="w")
        row += 1
        
        # Botón
        ttk.Button(scrollable_frame, text="🔬 CALCULAR", command=self.calculate).grid(row=row, column=0, columnspan=3, pady=10, sticky="ew")
        
        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        
        # Panel derecho (Outputs)
        right_frame = ttk.Frame(main_frame)
        right_frame.grid(row=0, column=1, sticky="nsew", padx=(5, 0))
        
        # Notebook para pestañas
        self.notebook = ttk.Notebook(right_frame)
        self.notebook.pack(fill=tk.BOTH, expand=True)
        
        # Pestaña Resultados
        self.results_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.results_frame, text="📋 Resultados")
        
        # Pestaña Gráficos
        self.plots_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.plots_frame, text="📊 Gráficos Orbitales")
        
        # Pestaña Perturbaciones
        self.perturb_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.perturb_frame, text="🌊 Análisis Perturbaciones")
        
        # Configurar grid weights
        main_frame.columnconfigure(1, weight=1)
        main_frame.rowconfigure(0, weight=1)
    
    def apply_preset(self):
        """Aplica configuraciones predefinidas de órbitas"""
        preset = self.orbit_preset.get()
        
        if preset == "LEO":
            self.alt_p_inicial.set(400.0)
            self.alt_a_inicial.set(400.0)
            self.alt_target.set(400.0)
            self.inclination_initial.set(51.6)  # ISS
            self.inclination_target.set(51.6)
            self.include_drag.set(True)
            
        elif preset == "MEO":
            self.alt_p_inicial.set(400.0)
            self.alt_a_inicial.set(400.0)
            self.alt_target.set(20200.0)  # GPS
            self.inclination_initial.set(28.5)
            self.inclination_target.set(55.0)  # GPS
            self.include_drag.set(False)
            
        elif preset == "GTO":
            self.alt_p_inicial.set(200.0)
            self.alt_a_inicial.set(35786.0)  # GTO típico
            self.alt_target.set(35786.0)
            self.inclination_initial.set(28.5)
            self.inclination_target.set(0.0)
            self.include_drag.set(False)
            
        elif preset == "GEO":
            self.alt_p_inicial.set(400.0)
            self.alt_a_inicial.set(400.0)
            self.alt_target.set(35786.0)
            self.inclination_initial.set(28.5)
            self.inclination_target.set(0.0)
            self.include_drag.set(False)
            self.include_sun.set(True)
            self.include_moon.set(True)
            
        elif preset == "LUNAR":
            self.alt_p_inicial.set(400.0)
            self.alt_a_inicial.set(400.0)
            self.alt_target.set(MOON_DIST - R_EARTH_KM)
            self.inclination_initial.set(28.5)
            self.inclination_target.set(28.5)
            self.include_drag.set(False)
            self.include_sun.set(True)
            self.include_moon.set(True)
        
    def calculate(self):
        try:
            # Validaciones
            if self.alt_p_inicial.get() > self.alt_a_inicial.get():
                messagebox.showerror("Error", "La altitud del perigeo no puede ser mayor que la del apogeo")
                return
            
            # Parámetros orbitales
            r_p_inicial = R_EARTH_KM + self.alt_p_inicial.get()
            r_a_inicial = R_EARTH_KM + self.alt_a_inicial.get()
            a_inicial = (r_p_inicial + r_a_inicial) / 2
            ecc_inicial = (r_a_inicial - r_p_inicial) / (r_a_inicial + r_p_inicial)
            period_inicial = calc_period(a_inicial)
            
            r_target = R_EARTH_KM + self.alt_target.get()
            period_target = calc_period(r_target)
            
            # Clasificación de órbitas
            orbit_class_initial = get_orbit_classification((self.alt_p_inicial.get() + self.alt_a_inicial.get()) / 2)
            orbit_class_target = get_orbit_classification(self.alt_target.get())
            
            # Calcular maniobra de transferencia orbital
            if self.maneuver_type.get() == "Hohmann":
                dv1, dv2, dv_total_transfer, T_transfer = hohmann_transfer(r_p_inicial, r_target)
                maneuver_details = f"Impulso 1 (Perigeo): {dv1*1000:.2f} m/s\nImpulso 2 (Apogeo): {dv2*1000:.2f} m/s"
                num_impulses = 2
            else:
                r_inter = R_EARTH_KM + self.alt_intermediate.get()
                dv1, dv2, dv3, dv_total_transfer, T_transfer = bielliptic_transfer(r_p_inicial, r_target, r_inter)
                maneuver_details = f"Impulso 1: {dv1*1000:.2f} m/s\nImpulso 2: {dv2*1000:.2f} m/s\nImpulso 3: {dv3*1000:.2f} m/s"
                num_impulses = 3
            
            # Calcular cambio de plano orbital si es necesario
            delta_i = abs(self.inclination_target.get() - self.inclination_initial.get())
            dv_plane_change = 0
            plane_change_details = ""
            
            if delta_i > 0.01:
                v_target = calc_velocity(r_target, r_target)
                dv_plane_change = calc_plane_change_dv(v_target, np.radians(delta_i))
                plane_change_details = f"\n\n🔄 CAMBIO DE PLANO ORBITAL:\nCambio de inclinación: {delta_i:.2f}°\nΔV adicional: {dv_plane_change*1000:.2f} m/s\n(Aplicado en nodo ascendente/descendente)"
            
            dv_total = dv_total_transfer + dv_plane_change
            
            # Perturbaciones
            perturbation_text = ""
            perturbation_data = {}
            
            if self.include_j2.get():
                i_rad = np.radians(self.inclination_initial.get())
                dOmega_dt, domega_dt = calc_perturbation_j2(a_inicial, ecc_inicial, i_rad)
                
                dOmega_day = np.degrees(dOmega_dt * 86400)
                domega_day = np.degrees(domega_dt * 86400)
                
                perturbation_data['j2'] = {
                    'dOmega_day': dOmega_day,
                    'domega_day': domega_day
                }
                
                perturbation_text += f"\n🔍 Perturbación J2 (Aplanamiento):\n"
                perturbation_text += f"  Precesión nodal: {dOmega_day:.4f} °/día\n"
                perturbation_text += f"  Precesión perigeo: {domega_day:.4f} °/día\n"
                perturbation_text += f"  (Acumulado 1 año): {dOmega_day*365:.2f}° (nodo), {domega_day*365:.2f}° (perigeo)\n"
                if self.alt_target.get() > 100000:
                    perturbation_text += f"  ⚠️ Nota: J2 tiene efecto reducido en órbitas muy altas\n"
            
            if self.include_drag.get():
                h_avg = (self.alt_p_inicial.get() + self.alt_a_inicial.get()) / 2
                if h_avg < 1000:  # Solo calcular si está en rango atmosférico
                    rho = atmospheric_density(h_avg)
                    v_avg = calc_velocity((r_p_inicial + r_a_inicial)/2, a_inicial)
                    drag_accel = 0.5 * rho * v_avg**2 * self.drag_coef.get() * self.cross_section.get() / self.sat_mass.get()
                    
                    energy_loss_per_orbit = drag_accel * v_avg * period_inicial / 1000
                    specific_energy = -MU_EARTH / (2 * a_inicial)
                    decay_time_years = abs(specific_energy / (energy_loss_per_orbit * (86400 * 365 / period_inicial)))
                    
                    perturbation_data['drag'] = {
                        'accel': drag_accel,
                        'energy_loss': energy_loss_per_orbit,
                        'lifetime': decay_time_years
                    }
                    
                    perturbation_text += f"\n💨 Arrastre Atmosférico:\n"
                    perturbation_text += f"  Densidad a {h_avg:.0f} km: {rho:.2e} kg/m³\n"
                    perturbation_text += f"  Aceleración: {drag_accel*1e6:.3f} μm/s²\n"
                    perturbation_text += f"  Pérdida energía/órbita: {energy_loss_per_orbit:.3e} km²/s²\n"
                    perturbation_text += f"  Vida útil estimada: {decay_time_years:.1f} años\n"
                else:
                    perturbation_text += f"\n💨 Arrastre Atmosférico:\n"
                    perturbation_text += f"  ⚠️ Órbita demasiado alta ({h_avg:.0f} km) - arrastre despreciable\n"
            
            if self.include_sun.get():
                sun_accel_max = MU_SUN * r_target / AU**3
                perturbation_data['sun'] = {'accel': sun_accel_max}
                perturbation_text += f"\n☀️ Perturbación Solar:\n"
                perturbation_text += f"  Aceleración máxima: {sun_accel_max*1e6:.3f} μm/s²\n"
                perturbation_text += f"  Efecto: Variación secular de excentricidad e inclinación\n"
            
            if self.include_moon.get():
                moon_accel_max = MU_MOON * r_target / MOON_DIST**3
                perturbation_data['moon'] = {'accel': moon_accel_max}
                perturbation_text += f"\n🌙 Perturbación Lunar:\n"
                perturbation_text += f"  Aceleración máxima: {moon_accel_max*1e6:.3f} μm/s²\n"
                perturbation_text += f"  Efecto: Variación periódica de elementos orbitales\n"
            
            # Combustible
            isp = 300
            g0 = 9.81 / 1000
            mass_ratio = np.exp(dv_total / (isp * g0))
            fuel_mass = self.sat_mass.get() * (mass_ratio - 1)
            
            # Guardar resultados
            self.current_results = {
                'maneuver_details': maneuver_details,
                'plane_change_details': plane_change_details,
                'dv_total': dv_total,
                'dv_transfer': dv_total_transfer,
                'dv_plane': dv_plane_change,
                'T_transfer': T_transfer,
                'fuel_mass': fuel_mass,
                'period_inicial': period_inicial,
                'period_target': period_target,
                'ecc_inicial': ecc_inicial,
                'a_inicial': a_inicial,
                'r_target': r_target,
                'perturbation_text': perturbation_text,
                'perturbation_data': perturbation_data,
                'num_impulses': num_impulses,
                'orbit_class_initial': orbit_class_initial,
                'orbit_class_target': orbit_class_target
            }
            
            # Mostrar resultados
            self.show_results()
            
            # Generar gráficos
            self.generate_plots(r_p_inicial, r_a_inicial, r_target, a_inicial, 
                              ecc_inicial, self.alt_target.get())
            
            # Generar análisis de perturbaciones
            self.generate_perturbation_analysis()
            
        except Exception as e:
            messagebox.showerror("Error", f"Error en el cálculo: {str(e)}")
    
    def show_results(self):
        if not self.current_results:
            return
            
        # Limpiar frame
        for widget in self.results_frame.winfo_children():
            widget.destroy()
        
        # Crear texto de resultados con scroll
        text_frame = ttk.Frame(self.results_frame)
        text_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        results_text = tk.Text(text_frame, wrap=tk.WORD, font=('Courier', 9))
        scrollbar = ttk.Scrollbar(text_frame, orient="vertical", command=results_text.yview)
        results_text.configure(yscrollcommand=scrollbar.set)
        
        results_text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        
        results = f"""
╔═══════════════════════════════════════════════════════════════╗
║        RESULTADOS DE LA CORRECCIÓN ORBITAL COMPLETA           ║
╚═══════════════════════════════════════════════════════════════╝

📡 CLASIFICACIÓN DE ÓRBITAS
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Órbita Inicial:  {self.current_results['orbit_class_initial']}
Órbita Objetivo: {self.current_results['orbit_class_target']}

🚀 DELTA-V REQUERIDO
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Tipo de Maniobra: {self.maneuver_type.get()}
{self.current_results['maneuver_details']}
{self.current_results['plane_change_details']}

🔥 Delta-V de Transferencia: {self.current_results['dv_transfer']*1000:.2f} m/s
🔄 Delta-V de Cambio de Plano: {self.current_results['dv_plane']*1000:.2f} m/s
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
💥 DELTA-V TOTAL: {self.current_results['dv_total']*1000:.2f} m/s
⏱️  Tiempo de transferencia: {self.current_results['T_transfer']/3600:.2f} horas ({self.current_results['T_transfer']/86400:.2f} días)

⛽ COMBUSTIBLE
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Masa del satélite: {self.sat_mass.get():.1f} kg
Masa de combustible: {self.current_results['fuel_mass']:.2f} kg
Masa total (sat + combustible): {self.sat_mass.get() + self.current_results['fuel_mass']:.2f} kg
Proporción combustible/satélite: {(self.current_results['fuel_mass']/self.sat_mass.get())*100:.1f}%
(Asumiendo ISP = 300s para propulsores químicos)

🪐 PARÁMETROS ORBITALES
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
INICIAL:
  Semieje mayor:    {self.current_results['a_inicial']:.1f} km
  Excentricidad:    {self.current_results['ecc_inicial']:.6f}
  Inclinación:      {self.inclination_initial.get():.2f}°
  Periodo:          {self.current_results['period_inicial']/3600:.2f} horas ({self.current_results['period_inicial']/60:.1f} min)
  Perigeo:          {self.alt_p_inicial.get():.1f} km
  Apogeo:           {self.alt_a_inicial.get():.1f} km

OBJETIVO:
  Semieje mayor:    {self.current_results['r_target']:.1f} km
  Excentricidad:    0.0000 (circular)
  Inclinación:      {self.inclination_target.get():.2f}°
  Periodo:          {self.current_results['period_target']/3600:.2f} horas ({self.current_results['period_target']/86400:.2f} días)
  Altitud:          {self.alt_target.get():.1f} km

🌍 PERTURBACIONES ORBITALES
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
{self.current_results['perturbation_text'] if self.current_results['perturbation_text'] else 'Ninguna perturbación incluida en el análisis'}

💡 RECOMENDACIONES
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
"""
        
        # Añadir recomendaciones basadas en los resultados
        if self.current_results['dv_total'] > 3.0:
            results += "⚠️  Alto consumo de Delta-V. Considere transferencia bi-elíptica.\n"
        
        if self.current_results['dv_plane'] > self.current_results['dv_transfer']:
            results += "⚠️  El cambio de plano domina el costo. Considere maniobra combinada.\n"
        
        if self.alt_p_inicial.get() < 300 and self.include_drag.get():
            results += "⚠️  Alta degradación por arrastre. Se requiere mantenimiento frecuente.\n"
        
        if self.current_results['fuel_mass']/self.sat_mass.get() > 0.5:
            results += "⚠️  Proporción de combustible muy alta. Verifique viabilidad de misión.\n"
        
        if self.alt_target.get() > 100000:
            results += "⚠️  Órbita muy alta. Perturbaciones lunares/solares son significativas.\n"
        
        if self.alt_target.get() > MOON_DIST * 0.5:
            results += "⚠️  Órbita cislunar. Considere modelo de 3 cuerpos (Tierra-Luna-Sat).\n"
        
        results += "\n✓ Cálculo completado exitosamente"
        
        results_text.insert('1.0', results)
        results_text.config(state='disabled')
        
    def generate_plots(self, r_p, r_a, r_target, a_inicial, ecc_inicial, alt_target):
        # Limpiar frame
        for widget in self.plots_frame.winfo_children():
            widget.destroy()
        
        # Crear figura con 4 subplots
        fig = Figure(figsize=(14, 10))
        
        # Gráfico 1: Vista orbital 2D con escala adaptativa
        ax1 = fig.add_subplot(221)
        
        # Órbita inicial
        nus = np.linspace(0, 2*np.pi, 200)
        rs_inicial = a_inicial * (1 - ecc_inicial**2) / (1 + ecc_inicial * np.cos(nus))
        xs_inicial = rs_inicial * np.cos(nus)
        ys_inicial = rs_inicial * np.sin(nus)
        ax1.plot(xs_inicial, ys_inicial, 'b-', label='Órbita Inicial', linewidth=2)
        
        # Tierra (primero para que quede atrás)
        circle_earth = plt.Circle((0, 0), R_EARTH_KM, color='lightblue', label='Tierra', zorder=5)
        ax1.add_patch(circle_earth)
        
        # Órbita objetivo
        circle_target = plt.Circle((0, 0), r_target, fill=False, color='g', linestyle='--', linewidth=2, label='Órbita Objetivo', zorder=6)
        ax1.add_patch(circle_target)
        
        # Órbita de transferencia y puntos ΔV
        if self.maneuver_type.get() == "Hohmann":
            a_transfer = (r_p + r_target) / 2
            ecc_transfer = (r_target - r_p) / (r_target + r_p)
            rs_transfer = a_transfer * (1 - ecc_transfer**2) / (1 + ecc_transfer * np.cos(nus))
            xs_transfer = rs_transfer * np.cos(nus)
            ys_transfer = rs_transfer * np.sin(nus)
            ax1.plot(xs_transfer, ys_transfer, 'r:', label='Transferencia Hohmann', linewidth=2, alpha=0.7, zorder=7)
            
            # Puntos de maniobra Hohmann (en posiciones correctas)
            # ΔV1 en el perigeo (ángulo 0)
            ax1.plot(r_p, 0, 'ro', markersize=12, label='ΔV1 (Perigeo)', zorder=15, markeredgecolor='darkred', markeredgewidth=2)
            # ΔV2 en el apogeo (ángulo π)
            ax1.plot(-r_target, 0, 'go', markersize=12, label='ΔV2 (Apogeo)', zorder=15, markeredgecolor='darkgreen', markeredgewidth=2)
            
        else:  # Bi-elíptica
            r_inter = R_EARTH_KM + self.alt_intermediate.get()
            
            # Primera elipse: r_p -> r_inter
            a_transfer1 = (r_p + r_inter) / 2
            ecc_transfer1 = (r_inter - r_p) / (r_inter + r_p)
            rs_transfer1 = a_transfer1 * (1 - ecc_transfer1**2) / (1 + ecc_transfer1 * np.cos(nus))
            xs_transfer1 = rs_transfer1 * np.cos(nus)
            ys_transfer1 = rs_transfer1 * np.sin(nus)
            ax1.plot(xs_transfer1, ys_transfer1, 'r:', label='Transferencia 1', linewidth=2, alpha=0.7, zorder=7)
            
            # Segunda elipse: r_inter -> r_target
            a_transfer2 = (r_inter + r_target) / 2
            ecc_transfer2 = abs(r_target - r_inter) / (r_target + r_inter)
            
            # Ajustar para que la segunda elipse conecte correctamente
            # El perigeo de la segunda elipse está en -r_inter
            rs_transfer2 = a_transfer2 * (1 - ecc_transfer2**2) / (1 + ecc_transfer2 * np.cos(nus))
            xs_transfer2 = -rs_transfer2 * np.cos(nus)
            ys_transfer2 = rs_transfer2 * np.sin(nus)
            ax1.plot(xs_transfer2, ys_transfer2, 'm:', label='Transferencia 2', linewidth=2, alpha=0.7, zorder=8)
            
            # Puntos de maniobra Bi-elíptica
            ax1.plot(r_p, 0, 'ro', markersize=12, label='ΔV1 (Perigeo)', zorder=15, markeredgecolor='darkred', markeredgewidth=2)
            ax1.plot(-r_inter, 0, 'yo', markersize=12, label='ΔV2 (Intermedia)', zorder=15, markeredgecolor='orange', markeredgewidth=2)
            ax1.plot(-r_target, 0, 'go', markersize=12, label='ΔV3 (Objetivo)', zorder=15, markeredgecolor='darkgreen', markeredgewidth=2)
            
            # Círculo de referencia intermedio
            circle_inter = plt.Circle((0, 0), r_inter, fill=False, color='yellow', linestyle=':', linewidth=1, alpha=0.5, label='Órbita Intermedia', zorder=6)
            ax1.add_patch(circle_inter)
        
        # Referencia órbita lunar si es relevante
        if alt_target > 100000:
            circle_moon_orbit = plt.Circle((0, 0), MOON_DIST, fill=False, color='gray', 
                                          linestyle=':', linewidth=1, alpha=0.5, label='Órbita Lunar')
            ax1.add_patch(circle_moon_orbit)
        
        # Referencia GEO si es relevante
        if 10000 < alt_target < 100000:
            circle_geo = plt.Circle((0, 0), R_EARTH_KM + GEO_ALT, fill=False, color='orange', 
                                   linestyle=':', linewidth=1, alpha=0.5, label='GEO')
            ax1.add_patch(circle_geo)
        
        ax1.set_xlabel('X (km)', fontsize=10)
        ax1.set_ylabel('Y (km)', fontsize=10)
        ax1.set_aspect('equal')
        ax1.grid(True, alpha=0.3)
        ax1.legend(loc='upper right', fontsize=8)
        ax1.set_title('Vista Orbital 2D', fontsize=12, fontweight='bold')
        
        # Gráfico 2: Altitud vs Anomalía
        ax2 = fig.add_subplot(222)
        
        anomalies_deg = np.linspace(0, 360, 400)
        anomalies_rad = np.radians(anomalies_deg)
        
        r_inicial_plot = a_inicial * (1 - ecc_inicial**2) / (1 + ecc_inicial * np.cos(anomalies_rad))
        alt_inicial_plot = r_inicial_plot - R_EARTH_KM
        
        ax2.plot(anomalies_deg, alt_inicial_plot, 'b-', label='Altitud Inicial', linewidth=2)
        ax2.axhline(y=alt_target, color='g', linestyle='--', linewidth=2, label=f'Objetivo ({alt_target:.0f} km)')
        ax2.fill_between(anomalies_deg, alt_inicial_plot, alt_target, alpha=0.2, color='orange')
        
        # Referencias orbitales
        if alt_target > 10000:
            ax2.axhline(y=GEO_ALT, color='orange', linestyle=':', linewidth=1, alpha=0.5, label='GEO')
        
        ax2.set_xlabel('Anomalía Verdadera (grados)', fontsize=10)
        ax2.set_ylabel('Altitud (km)', fontsize=10)
        ax2.legend(fontsize=8)
        ax2.grid(True, alpha=0.3)
        ax2.set_title('Variación de Altitud', fontsize=12, fontweight='bold')
        
        # Gráfico 3: Velocidad orbital
        ax3 = fig.add_subplot(223)
        
        v_inicial = np.sqrt(MU_EARTH * (2/r_inicial_plot - 1/a_inicial))
        v_target = np.sqrt(MU_EARTH / r_target)
        
        ax3.plot(anomalies_deg, v_inicial, 'b-', label='Velocidad Inicial', linewidth=2)
        ax3.axhline(y=v_target, color='g', linestyle='--', linewidth=2, label=f'Velocidad Objetivo')
        
        ax3.set_xlabel('Anomalía Verdadera (grados)', fontsize=10)
        ax3.set_ylabel('Velocidad (km/s)', fontsize=10)
        ax3.legend(fontsize=8)
        ax3.grid(True, alpha=0.3)
        ax3.set_title('Perfil de Velocidad Orbital', fontsize=12, fontweight='bold')
        
        # Gráfico 4: Energía específica
        ax4 = fig.add_subplot(224)
        
        energy_inicial = v_inicial**2 / 2 - MU_EARTH / r_inicial_plot
        energy_target = v_target**2 / 2 - MU_EARTH / r_target
        
        ax4.plot(anomalies_deg, energy_inicial, 'b-', label='Energía Inicial', linewidth=2)
        ax4.axhline(y=energy_target, color='g', linestyle='--', linewidth=2, label=f'Energía Objetivo')
        
        ax4.set_xlabel('Anomalía Verdadera (grados)', fontsize=10)
        ax4.set_ylabel('Energía Específica (km²/s²)', fontsize=10)
        ax4.legend(fontsize=8)
        ax4.grid(True, alpha=0.3)
        ax4.set_title('Energía Orbital', fontsize=12, fontweight='bold')
        
        fig.tight_layout()
        
        # Incrustar en tkinter
        canvas = FigureCanvasTkAgg(fig, master=self.plots_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
    def generate_perturbation_analysis(self):
        # Limpiar frame
        for widget in self.perturb_frame.winfo_children():
            widget.destroy()
        
        if not self.current_results or not self.current_results['perturbation_data']:
            ttk.Label(self.perturb_frame, text="No hay datos de perturbaciones disponibles.\nActive al menos una perturbación y calcule.", 
                     font=('Arial', 12)).pack(pady=50)
            return
        
        # Crear figura
        fig = Figure(figsize=(14, 10))
        
        perturb_data = self.current_results['perturbation_data']
        
        # Gráfico de efectos acumulativos
        days_display = 30 if self.alt_target.get() < 10000 else 365
        
        if 'j2' in perturb_data:
            ax1 = fig.add_subplot(221)
            days = np.linspace(0, days_display, 100)
            
            dOmega_accum = perturb_data['j2']['dOmega_day'] * days
            domega_accum = perturb_data['j2']['domega_day'] * days
            
            ax1.plot(days, dOmega_accum, 'b-', label='Nodo Ascendente (Ω)', linewidth=2)
            ax1.plot(days, domega_accum, 'r-', label='Arg. Perigeo (ω)', linewidth=2)
            
            ax1.set_xlabel('Tiempo (días)', fontsize=10)
            ax1.set_ylabel('Precesión Acumulada (°)', fontsize=10)
            ax1.legend(fontsize=9)
            ax1.grid(True, alpha=0.3)
            ax1.set_title('Efecto J2: Precesión Acumulada', fontsize=11, fontweight='bold')
        
        # Gráfico de degradación orbital por arrastre
        if 'drag' in perturb_data:
            ax2 = fig.add_subplot(222)
            days = np.linspace(0, min(days_display, perturb_data['drag']['lifetime']*365), 100)
            
            a_inicial = self.current_results['a_inicial']
            energy_loss_per_day = perturb_data['drag']['energy_loss'] * (86400 / self.current_results['period_inicial'])
            
            a_decay = []
            for d in days:
                delta_E = energy_loss_per_day * d
                a_new = -MU_EARTH / (2 * (-MU_EARTH/(2*a_inicial) - delta_E))
                a_decay.append(a_new)
            
            ax2.plot(days, np.array(a_decay) - R_EARTH_KM, 'orange', linewidth=2)
            ax2.axhline(y=100, color='r', linestyle='--', label='Límite reentrada (~100 km)', linewidth=1)
            
            ax2.set_xlabel('Tiempo (días)', fontsize=10)
            ax2.set_ylabel('Altitud del Semieje Mayor (km)', fontsize=10)
            ax2.legend(fontsize=9)
            ax2.grid(True, alpha=0.3)
            ax2.set_title('Degradación Orbital por Arrastre', fontsize=11, fontweight='bold')
        
        # Gráfico comparativo de magnitudes de perturbaciones
        ax3 = fig.add_subplot(223)
        
        perturb_names = []
        perturb_accels = []
        
        if 'drag' in perturb_data:
            perturb_names.append('Arrastre')
            perturb_accels.append(perturb_data['drag']['accel'] * 1e6)
        
        if 'sun' in perturb_data:
            perturb_names.append('Sol')
            perturb_accels.append(perturb_data['sun']['accel'] * 1e6)
        
        if 'moon' in perturb_data:
            perturb_names.append('Luna')
            perturb_accels.append(perturb_data['moon']['accel'] * 1e6)
        
        if 'j2' in perturb_data:
            j2_effect = abs(perturb_data['j2']['dOmega_day']) * 1000
            perturb_names.append('J2')
            perturb_accels.append(j2_effect)
        
        if perturb_names:
            colors = ['orange', 'gold', 'silver', 'skyblue'][:len(perturb_names)]
            bars = ax3.barh(perturb_names, perturb_accels, color=colors)
            
            ax3.set_xlabel('Magnitud del Efecto (μm/s² o equiv.)', fontsize=10)
            ax3.set_title('Comparación de Perturbaciones', fontsize=11, fontweight='bold')
            ax3.grid(True, alpha=0.3, axis='x')
            
            for i, (bar, val) in enumerate(zip(bars, perturb_accels)):
                ax3.text(val, i, f' {val:.3f}', va='center', fontsize=9)
        
        # Tabla resumen
        ax4 = fig.add_subplot(224)
        ax4.axis('off')
        
        table_data = []
        table_data.append(['Perturbación', 'Efecto Principal', 'Magnitud'])
        table_data.append(['─' * 13, '─' * 25, '─' * 16])
        
        if 'j2' in perturb_data:
            table_data.append(['J2', 'Precesión nodal', f"{perturb_data['j2']['dOmega_day']:.4f} °/día"])
            table_data.append(['', 'Precesión perigeo', f"{perturb_data['j2']['domega_day']:.4f} °/día"])
        
        if 'drag' in perturb_data:
            table_data.append(['Arrastre', 'Decaimiento orbital', f"{perturb_data['drag']['lifetime']:.1f} años"])
        
        if 'sun' in perturb_data:
            table_data.append(['Solar', 'Variación secular e, i', f"{perturb_data['sun']['accel']*1e6:.3f} μm/s²"])
        
        if 'moon' in perturb_data:
            table_data.append(['Lunar', 'Variación periódica', f"{perturb_data['moon']['accel']*1e6:.3f} μm/s²"])
        
        table = ax4.table(cellText=table_data, cellLoc='left', loc='center',
                         colWidths=[0.25, 0.45, 0.3])
        table.auto_set_font_size(False)
        table.set_fontsize(9)
        table.scale(1, 2)
        
        for i in range(len(table_data)):
            for j in range(3):
                cell = table[(i, j)]
                if i == 0:
                    cell.set_facecolor('#4CAF50')
                    cell.set_text_props(weight='bold', color='white')
                elif i == 1:
                    cell.set_facecolor('#f0f0f0')
                else:
                    cell.set_facecolor('white' if i % 2 == 0 else '#f9f9f9')
        
        ax4.set_title('Resumen de Perturbaciones', fontsize=11, fontweight='bold', pad=20)
        
        fig.tight_layout()
        
        canvas = FigureCanvasTkAgg(fig, master=self.perturb_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

if __name__ == "__main__":
    root = tk.Tk()
    app = OrbitalCalculatorApp(root)
    root.mainloop()
