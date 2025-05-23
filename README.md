# 🧭 Polygon Orientation Analyzer

This project analyzes building footprints or land parcel polygons from a shapefile and calculates the orientation of their exterior sides that do not touch neighboring polygons. It identifies predominant facade directions such as **Norte (North)**, **Este (East)**, **Sur (South)**, etc., and outputs results to CSV files.

---

## 📂 Features

- Calculates **angles and cardinal orientations** of exterior polygon sides.
- Filters out shared or touching boundaries with adjacent geometries.
- Supports **MultiPolygon** and **Polygon** geometries.
- Outputs:
  - Orientation summary per polygon
  - Details of each detected exterior side (coordinates, angle, direction, length)

---

## 📦 Requirements

Install dependencies using:

```bash
pip install -r requirements.txt
