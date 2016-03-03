from curve_augmentation import get_transformations

d = {

    'rotation': {

        'theta': [3, 5, 7, 10],
        'center': [-1]
    },
    'noise': {

        'sigma': [0.5, 1., 1.5, 2.],
        'seed': [1]
    }
}

print get_transformations(d, 16)