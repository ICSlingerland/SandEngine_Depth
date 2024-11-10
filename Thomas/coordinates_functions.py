import numpy as np
import pyproj
import matplotlib.pyplot as plt
def convert_rd_to_wgs84(points_rd):
    # Define the RD New (Amersfoort) CRS
    proj_rd_new = pyproj.Proj(init="epsg:28992")

    # Define the WGS84 CRS (used for lat/lon)
    proj_wgs84 = pyproj.Proj(init="epsg:4326")

    # Convert RD New to WGS84 (Lat, Lon) and retain height
    points_wgs84 = [
        (*pyproj.transform(proj_rd_new, proj_wgs84, point[0], point[1]), point[2])
        for point in points_rd
    ]
    
    return points_wgs84

from PIL import Image

def get_image_coordinates(img_path, center_x, center_y, block_size=10, color=(255, 0, 0)):
    """
    Highlights a square area of specified size around a center pixel in the given color.
    
    :param img_path: Path to the image file
    :param center_x: X coordinate of the center pixel
    :param center_y: Y coordinate of the center pixel
    :param block_size: Size of the square block to color (default is 10x10)
    :param color: RGB color to fill the block (default is red)
    """
    # Open the image
    img = Image.open(img_path)
    
    # Calculate top-left corner of the block
    pix_x = center_x - block_size // 2
    pix_y = center_y - block_size // 2
    
    # Modify a block_size x block_size area around the center pixel
    for i in range(pix_x, pix_x + block_size):
        for j in range(pix_y, pix_y + block_size):
            img.putpixel((i, j), color)
    
    # Show the modified image
    img.show()
    # Optionally save or display the image with matplotlib
    # plt.imshow(img)

# Example usage




def compute_projection_matrix(points_rd_array, image_points):
    """
    Computes the projection matrix P using real-world coordinates and image coordinates.
    
    :param points_rd_array: Array of real-world coordinates, shape (N, 3) where N is the number of points
    :param image_points: Array of image coordinates, shape (N, 2)
    :return: Projection matrix P of shape (3, 4)
    """
    # Number of points
    N = points_rd_array.shape[0]
    
    # Initialize A matrix and b vector
    A = []
    b = []
    
    # Construct the A matrix and b vector based on the provided points
    for i in range(N):
        X, Y, Z = points_rd_array[i, :]
        x_prime, y_prime = image_points[i, :]
        
        # Append rows for the A matrix
        A.append([X, Y, Z, 1, 0, 0, 0, 0, 0, 0, 0, 0])
        A.append([0, 0, 0, 0, X, Y, Z, 1, 0, 0, 0, 0])
        A.append([0, 0, 0, 0, 0, 0, 0, 0, X, Y, Z, 1])
        
        # Append corresponding values for the b vector
        b.append(x_prime)
        b.append(y_prime)
        b.append(1)
    
    # Convert A and b to numpy arrays
    A = np.array(A)
    b = np.array(b)
    
    # Solve the least-squares problem to find the projection parameters p
    p, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
    
    # Reshape the solution into a 3x4 matrix
    P = p.reshape(3, 4)
    
    return P

def plot_gcp_residuals(gcp_array, residuals, padding_x, padding_y, title="GCP Residuals Plot", x_label="Easting direction", y_label="Northing direction", scale=2, resolution=(8, 6)):
    """
    Plots the Ground Control Points (GCP) residuals along with the Argus Tower position on a quiver plot.
    
    :param gcp_array: Array of GCP coordinates (easting, northing, elevation), shape (N, 3)
    :param residuals: Array of residual errors (easting_error, northing_error), shape (N, 2)
    :param argus_coordinates: Tuple with the coordinates of the Argus Tower (easting, northing)
    :param padding_x: Padding added to x-axis limits
    :param padding_y: Padding added to y-axis limits
    :param title: Title of the plot
    :param x_label: Label for the x-axis
    :param y_label: Label for the y-axis
    :param scale: Scale for the residuals in the quiver plot (default is 15)
    :param resolution: Resolution of the plot figure in inches (width, height), default (8, 6)
    :return: Displays the plot and outputs the image
    """
    # Extract easting and northing components for GCPs
    easting = gcp_array[:, 0]
    northing = gcp_array[:, 1]

    # Extract residual errors for easting and northing
    easting_errors = residuals[:, 0]
    northing_errors = residuals[:, 1]

    # Define GCP labels
    gcp_labels = ["GCP4", "GCP5", "GCP6", "GCP7", "GCP8", "GCP9", "GCP10"]

    # Calculate axis limits with padding
    min_easting, max_easting = np.min(easting) - padding_x, np.max(easting) + padding_x
    min_northing, max_northing = np.min(northing) - padding_y, np.max(northing) + padding_y

    # Define offset for the labels
    x_offset = 6
    y_offset = 6



    # Create a quiver plot
    plt.figure(figsize=resolution)
    plt.quiver(easting, northing, easting_errors, northing_errors, angles='xy', scale_units='xy', scale=scale, color='blue')
    plt.scatter(easting, northing, color='red', label='GCPs')  # Mark GCP locations

    # Set axis limits with padding
    plt.xlim(min_easting, max_easting)
    plt.ylim(min_northing, max_northing)

    # Flip the y-axis and x-axis if needed
    plt.gca().invert_yaxis()
   

    # Add labels for each GCP with offset
    for i, label in enumerate(gcp_labels):
        plt.text(easting[i] + x_offset, northing[i] + y_offset, label, fontsize=10, ha='right', va='bottom')

    # Add labels and grid
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(title)
    plt.legend()
    plt.grid(True)

    # Display the plot
    plt.show()


def create_coordinate_grid(image_coordinates, lidar_points, resolution=10):
    """
    Create a grid over the given coordinates with specified resolution,
    and assign points to grid cells.

    Parameters:
    - image_coordinates (ndarray): An array of shape (N, 2) with x (easting) and y (northing) coordinates.
    - lidar_points (ndarray): An array of shape (N, 2) with real-world coordinates corresponding to image coordinates.
    - resolution (int): The resolution of the grid in pixels (default is 10).

    Returns:
    - dict: A dictionary with grid cells as keys (i, j) and values as dictionaries containing image and lidar points.
    """

    # Define grid boundaries
    min_E, max_E = np.min(image_coordinates[:, 0]), np.max(image_coordinates[:, 0])
    min_N, max_N = np.min(image_coordinates[:, 1]), np.max(image_coordinates[:, 1])
    
    # Create bins for the grid
    easting_bins = np.arange(min_E, max_E + resolution, resolution)
    northing_bins = np.arange(min_N, max_N + resolution, resolution)

    # Make the grid using meshgrid
    easting_grid, northing_grid = np.meshgrid(easting_bins, northing_bins)

    # Initialize a dictionary to store the points in each grid box
    grid = {}

    # Iterate over each grid cell
    for i in range(easting_grid.shape[0]):
        print(f"Processing row {i + 1}/{easting_grid.shape[0]}...")  # Display progress for each row
        for j in range(northing_grid.shape[1]):
            # Center and bounds of the current grid cell
            east_center = easting_grid[i, j]
            north_center = northing_grid[i, j]
            east_min_region = east_center - resolution / 2
            east_max_region = east_center + resolution / 2
            north_min_region = north_center - resolution / 2
            north_max_region = north_center + resolution / 2

            # Find indices of image coordinates within the current grid box
            image_coordinates_indices = np.where(
                (image_coordinates[:, 0] >= east_min_region) & (image_coordinates[:, 0] < east_max_region) &
                (image_coordinates[:, 1] >= north_min_region) & (image_coordinates[:, 1] < north_max_region)
            )[0]

            # Assign image and real coordinates to the current grid cell
            grid_key = (i, j)
            grid[grid_key] = {
                'image_coordinates': image_coordinates[image_coordinates_indices],
                'real_coordinates': lidar_points[image_coordinates_indices]
            }

    return grid

import numpy as np

def create_distance_grid(image_coordinates, lidar_points, argus_tower, resolution=5):
    """
    Creates a grid from lidar points and image coordinates, then calculates
    the minimum distance from a reference point (e.g., Argus Tower) for each cell.

    Parameters:
    - image_coordinates (ndarray): Array of shape (N, 2) with x (easting) and y (northing) coordinates.
    - lidar_points (ndarray): Array of shape (N, 2) with real-world coordinates corresponding to image coordinates.
    - argus_tower (ndarray): Array with the reference point coordinates to calculate distances.
    - resolution (int): The resolution of the grid in pixels (default is 5).

    Returns:
    - dict: A dictionary where each key is a grid cell, with values containing min distances to the reference point.
    """

    # Define grid boundaries
    min_E, max_E = np.min(image_coordinates[:, 0]), np.max(image_coordinates[:, 0])
    min_N, max_N = np.min(image_coordinates[:, 1]), np.max(image_coordinates[:, 1])

    # Create bins for the grid
    easting_bins = np.arange(min_E, max_E + resolution, resolution)
    northing_bins = np.arange(min_N, max_N + resolution, resolution)

    # Digitize to find which grid cell each point belongs to
    easting_indices = np.digitize(image_coordinates[:, 0], easting_bins) - 1
    northing_indices = np.digitize(image_coordinates[:, 1], northing_bins) - 1

    # Initialize the grid dictionary to store points and min distances
    grid = {}
    min_distances = {}

    # Assign each point to its corresponding grid cell
    for idx in range(len(image_coordinates)):
        i, j = northing_indices[idx], easting_indices[idx]  # Grid cell indices

        # Skip points out of the grid bounds
        if i < 0 or i >= len(northing_bins) - 1 or j < 0 or j >= len(easting_bins) - 1:
            continue

        # Create a unique key for the grid cell
        grid_key = (i, j)
        
        # Initialize the cell if it doesn't exist
        if grid_key not in grid:
            grid[grid_key] = {'image_coordinates': [], 'real_coordinates': []}
        
        # Append image and real coordinates to the grid cell
        grid[grid_key]['image_coordinates'].append(image_coordinates[idx])
        grid[grid_key]['real_coordinates'].append(lidar_points[idx])

    # Calculate minimum distance for each grid cell
    for key, value in grid.items():
        if value['real_coordinates']:  # Only process non-empty cells
            # Convert lists to arrays
            value['image_coordinates'] = np.array(value['image_coordinates'])
            value['real_coordinates'] = np.array(value['real_coordinates'])

            # Calculate distances to the reference point for all coordinates in the cell
            distances = np.linalg.norm(value['real_coordinates'] - argus_tower[:2], axis=1)

            # Find the minimum distance for the grid cell
            min_distances[key] = np.min(distances)

    print("Grid processing and distance calculation complete.")
    return grid, min_distances

