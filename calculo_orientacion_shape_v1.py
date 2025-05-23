import geopandas as gpd
from shapely.geometry import Polygon, MultiPolygon, LineString, Point
from shapely.affinity import translate
import math
import pandas as pd
import warnings
import string

# Define a buffer distance in meters 
BUFFER_DISTANCE = 11  # You can adjust this as needed

warnings.filterwarnings("ignore")

# Function to calculate the angle of a side in degrees with respect to cardinal points
def calcular_angulo(punto1, punto2):
    dx = punto2[0] - punto1[0]
    dy = punto2[1] - punto1[1]
    angle_rad = math.atan2(dy, dx)
    angle_deg = math.degrees(angle_rad)
    compass_angle = (90 - angle_deg) % 360
    return compass_angle

# Function to assign cardinal point based on the angle
def asignar_punto_cardinal(angle):
    directions = [
        ('Norte', 60, 120),
        ('Noroeste', 30, 60),
        ('Noreste', 120, 150),
        ('Sur', 240, 300),
        ('Suroeste', 300, 330),
        ('Sureste', 210, 240),
        ('Este', 150, 210),
        ('Oeste', 330, 360),
        ('Oeste', 0, 30),
    ]
    for direction, min_angle, max_angle in directions:
        if min_angle <= angle < max_angle:
            return direction
    return 'Indefinido'

# Función para calcular el vector perpendicular
def perpendicular_vector(line, distance):
    """
    Calcula un vector perpendicular al segmento de línea con una distancia específica.
    """
    x1, y1 = line.coords[0]
    x2, y2 = line.coords[1]
    
    # Dirección del lado (vector)
    dx = x2 - x1
    dy = y2 - y1
    
    # Longitud del vector
    length = math.sqrt(dx**2 + dy**2)
    
    # Vector unitario perpendicular (girado 90 grados)
    ux = -dy / length
    uy = dx / length
    
    # Vector perpendicular escalado a la distancia deseada
    px = ux * distance
    py = uy * distance
    
    return px, py

def procesar_lados_exteriores(polygon, index, parte, record_id, resultados, lados, otros_poligonos, gdfcrs):
    coords = list(polygon.exterior.coords)
    num_lados = len(coords) - 1  # Número de lados
    sin_sum = 0
    cos_sum = 0
    total_longitud = 0
    BUFFER_DISTANCE = 10  # Buffer de 10 metros

    print(f"Número de lados: {num_lados}")
    
    for i in range(num_lados):
        punto1 = coords[i]
        punto2 = coords[i + 1]
        
        print(f"Procesando lado {i + 1}: Punto1 = {punto1}, Punto2 = {punto2}")
        
        lado = LineString([punto1, punto2])
        
        # Crear un buffer lineal en lugar de circular alrededor del lado
        buffer_lado = lado.buffer(BUFFER_DISTANCE, cap_style=2)  # cap_style=2 asegura un buffer recto en los extremos

        # Verificar si los dos vértices tocan el otro polígono
        punto1_geom = Point(punto1)
        punto2_geom = Point(punto2)
        toca_punto1 = otros_poligonos.touches(punto1_geom)
        toca_punto2 = otros_poligonos.touches(punto2_geom)

        # Revisar si el lado cruza, toca o si el buffer interseca otros polígonos
        cruza_otro_poligono = otros_poligonos.crosses(lado)
        interseccion_buffer = otros_poligonos.intersects(buffer_lado)
        interseccion = otros_poligonos.intersects(lado)

        print(f"Lado {i+1}: Cruza={cruza_otro_poligono.any()}, Toca Puntos={toca_punto1.any() and toca_punto2.any()}, Intersección ={interseccion.any()}, Intersección Buffer={interseccion_buffer.any()}")

        # Excluir si ambos vértices tocan el polígono o si cruza/intersecta el buffer o el lado
        if (cruza_otro_poligono.any() == True) or (interseccion_buffer.any() == True) or (toca_punto1.any() == True and toca_punto2.any() == True):
            print(f"Lado {i+1} excluido por intersección o contacto.")
            continue

        # Calcular el ángulo del lado
        angle = calcular_angulo(punto1, punto2)
        direccion = asignar_punto_cardinal(angle)
        
        # Calcular la longitud del lado
        longitud = math.sqrt((punto2[0] - punto1[0])**2 + (punto2[1] - punto1[1])**2)
        total_longitud += longitud
        
        sin_sum += math.sin(math.radians(angle)) * longitud
        cos_sum += math.cos(math.radians(angle)) * longitud

        # Convertir punto1 a coordenadas geográficas (EPSG:4326)
        punto1_geom_4326 = gpd.GeoSeries([Point(punto1)], crs=gdfcrs).to_crs(epsg=4326).iloc[0]
        punto2_geom_4326 = gpd.GeoSeries([Point(punto2)], crs=gdfcrs).to_crs(epsg=4326).iloc[0]
        
        # Almacenar información del lado en la lista de resultados
        lados.append({
            'Polígono': record_id,
            'Parte': parte,
            'Lado': i + 1,
            'Lat_1': punto1_geom_4326.y,
            'Lon_1': punto1_geom_4326.x,
            'Lat_2': punto2_geom_4326.y,
            'Lon_2': punto2_geom_4326.x,
            'Ángulo (°)': round(angle, 2),
            'Orientación': direccion,
            'Longitud': round(longitud, 2)
        })
    
    # Procesar el último lado que cierra el polígono
    punto1 = coords[-1]
    punto2 = coords[0]
    print(f"Procesando lado final: Punto1 = {punto1}, Punto2 = {punto2}")

    
    # Calcular el ángulo promedio si hay longitud total
    if total_longitud > 0:
        angle_promedio_rad = math.atan2(sin_sum / total_longitud, cos_sum / total_longitud)
        angle_promedio_deg = math.degrees(angle_promedio_rad)
        compass_angle = (90 - angle_promedio_deg) % 360
        angle_promedio = compass_angle
    else:
        angle_promedio = 0

    # Diccionario para almacenar las longitudes por orientación
    orientaciones_longitud = {
        'Norte': 0,
        'Noreste': 0,
        'Este': 0,
        'Sureste': 0,
        'Sur': 0,
        'Suroeste': 0,
        'Oeste': 0,
        'Noroeste': 0
    }

    # Sumar las longitudes por orientación
    for lado in lados:
        if lado['Polígono'] == record_id:
            orientacion = lado['Orientación']
            longitud = lado['Longitud']
            orientaciones_longitud[orientacion] += longitud
    
    # Determinar la orientación predominante
    orientacion_promedio = max(orientaciones_longitud, key=orientaciones_longitud.get)        
    
    # Almacenar resultados en la lista final
    resultados.append({
        'Polígono': record_id,
        'Ángulo Promedio (°)': round(angle_promedio, 2),
        'Orientación': orientacion_promedio
    })

# Function to process geometry (Polygon or MultiPolygon)
def procesar_geometria_exterior(geometry, geometria_union, index, record_id, resultados, lados, otros_poligonos, gdfcrs):
    exterior = geometry.difference(geometria_union)
    
    if exterior.is_empty:
        print(f"Polígono {record_id}: No tiene fronteras exteriores únicas.")
        return

    if isinstance(exterior, MultiPolygon):
        for i, poly in enumerate(exterior.geoms):
            procesar_lados_exteriores(poly, index, f"{i + 1}", record_id, resultados, lados, otros_poligonos, gdfcrs)
    elif isinstance(exterior, Polygon):
        procesar_lados_exteriores(exterior, index, "1", record_id, resultados, lados, otros_poligonos, gdfcrs)
    else:
        print(f"Geometría no soportada en el índice {record_id}: {type(exterior)}")

# Main function to load the shapefile and process polygons based on references from CSV
def analizar_orientacion_exterior_shapefile(ruta_shapefile, ruta_csv_referencias, salida_orientacion_csv=None, salida_lados_csv=None):
    gdf = gpd.read_file(ruta_shapefile)

    # Convert CRS to UTM for consistent spatial operations
    utm_crs = gdf.estimate_utm_crs()
    gdf = gdf.to_crs(utm_crs)

    # Load references from the CSV
    referencias_df = pd.read_csv(ruta_csv_referencias, sep=";")
    referencias_df = referencias_df[['referencia_catastral', 'planta']].drop_duplicates()

    if gdf.empty or 'geometry' not in gdf.columns:
        raise ValueError("El archivo shapefile no contiene geometría válida.")

    resultados = []
    lados = []

    # Normalize the 'referencia' column in the GeoDataFrame
    gdf['referencia'] = gdf['referencia'].str.strip().str.upper()
    gdf['planta'] = gdf['planta'].astype(str).str.strip()

    # Iterate over all rows in the GeoDataFrame
    for index, row in gdf.iterrows():
        geometry = row.geometry
        record_id = row['referencia'].upper()
        planta = row['planta']
        print(f"Procesando Polígono {record_id} en planta {planta}...")

        # Create a GeoSeries of all other polygons
        otros_poligonos_planta = gdf.drop(gdf.index[index])
        otros_poligonos = otros_poligonos_planta[(otros_poligonos_planta['planta'] == planta)]
        if otros_poligonos.empty:
            print(f"Vacio.")
        # Union all other polygons into a single geometry
        union_otros = otros_poligonos.unary_union

        # Process geometry to obtain the outer boundaries that do not intersect
        procesar_geometria_exterior(geometry, union_otros, index, record_id, resultados, lados, otros_poligonos, gdf.crs)

    if resultados:
        # Convert results to a DataFrame
        df_resultados = pd.DataFrame(resultados)
        df_lados = pd.DataFrame(lados)

        # Save to CSV
        if salida_orientacion_csv:
            df_resultados.to_csv(salida_orientacion_csv, index=False, sep=";")
        if salida_lados_csv:
            df_lados.to_csv(salida_lados_csv, index=False, sep=";")
        print(df_resultados)
    else:
        print("No se encontraron lados exteriores únicos.")


# Path to the shapefile (.shp)
ruta_shapefile = "C:/Users/TechAir/Documents/GitHub/refcat_orientation_exterior/refcat_orientation_exterior/orientacion/sophiq_kml_con_rf/sophiq_kml_con_rf.shp"
# Path to the CSV file containing the references
ruta_csv_referencias = "C:/Users/TechAir/Documents/GitHub/refcat_orientation_exterior/refcat_orientation_exterior/orientacion/calculo_orientacion_shape_v1/referencias.csv"

# Optional: Save results to a CSV file
ruta_salida_orientacion_csv = 'C:/Users/TechAir/Documents/GitHub/refcat_orientation_exterior/refcat_orientation_exterior/orientacion/calculo_orientacion_shape_v1/orientacion_poligonos.csv' 
ruta_salida_lados_csv = 'C:/Users/TechAir/Documents/GitHub/refcat_orientation_exterior/refcat_orientation_exterior/orientacion/calculo_orientacion_shape_v1/orientacion_lados.csv'

# Execute the analysis and save the results
analizar_orientacion_exterior_shapefile(
    ruta_shapefile=ruta_shapefile,
    ruta_csv_referencias=ruta_csv_referencias,
    salida_orientacion_csv=ruta_salida_orientacion_csv,
    salida_lados_csv=ruta_salida_lados_csv
)
